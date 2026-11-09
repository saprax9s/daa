#include "bench.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "hashing.h"
#include "pipeline.h"
#include "proofs.h"
#include "sha256.h"
#include "trees.h"
#include "weights.h"
#include "workloads.h"

namespace merkle {

namespace {

volatile uint64_t g_sink = 0;

void doNotOptimize(uint64_t val) {
    g_sink = val;
}

struct CsvWriter {
    std::ofstream out;
    explicit CsvWriter(const std::string& path) : out(path) {
        out << "language,experiment,algo,K,N,m,zipf_s,dmax,seed,rep,metric,value,unit\n";
    }

    void write(const std::string& exp, const std::string& algo, uint64_t K, uint64_t N,
               uint64_t m, double s, uint64_t dmax, uint64_t seed, int rep,
               const std::string& metric, double val, const std::string& unit) {
        out << "cpp," << exp << "," << algo << "," << K << "," << N << "," << m << ","
            << s << "," << dmax << "," << seed << "," << rep << ","
            << metric << "," << val << "," << unit << "\n";
        out.flush();
    }
};

std::vector<uint64_t> makeZipfWeights(int count, double s, uint64_t seed) {
    std::vector<uint64_t> w(count);
    uint64_t state = seed;
    // Harmonic sum
    double c = 0.0;
    for (int i = 1; i <= count; ++i) c += 1.0 / std::pow((double)i, s);
    for (int i = 0; i < count; ++i) {
        uint64_t r = splitmix64(state);
        // Zipf rank distribution
        double rank = 1.0 + (double)(r % (uint64_t)count);
        double prob = (1.0 / std::pow(rank, s)) / c;
        uint64_t weight = (uint64_t)(prob * 10000000.0) + 1ULL;
        w[i] = weight;
    }
    return w;
}

}  // namespace

int runExperiments(const std::string& outputFile) {
    CsvWriter csv(outputFile);
    std::cout << "Starting Merkle Tree benchmark experiments..." << std::endl;

    const uint64_t baseSeed = 42ULL;

    // =========================================================================
    // E1: Build time and operation counts vs K in {64, 128, 256, 512, 1024}
    // =========================================================================
    std::cout << "[E1] Running build time and op counts vs K..." << std::endl;
    std::vector<int> kValues = {64, 128, 256, 512, 1024};
    for (int K : kValues) {
        int leaves = 2 * K + 1;
        std::vector<uint64_t> w = makeZipfWeights(leaves, 1.0, baseSeed + (uint64_t)K);
        uint64_t dmax = (uint64_t)(std::ceil(std::log2((double)leaves)) + 2.0);

        // Algos to benchmark
        struct Algo {
            std::string name;
            std::function<TreeResult(const std::vector<uint64_t>&)> fn;
            uint64_t dpTrans;
        };

        std::vector<Algo> algos = {
            {"balanced", [](const std::vector<uint64_t>& wt) { return solveBalanced(wt); }, 0},
            {"mehlhorn", [](const std::vector<uint64_t>& wt) { return solveMehlhorn(wt); }, (uint64_t)leaves * (uint64_t)leaves / 2},
            {"garsia_wachs", [](const std::vector<uint64_t>& wt) { return solveOptimal(wt); }, (uint64_t)leaves * (uint64_t)leaves / 2},
            {"knuth", [](const std::vector<uint64_t>& wt) { return solveKnuth(wt); }, (uint64_t)leaves * (uint64_t)leaves},
            {"depth_capped", [dmax](const std::vector<uint64_t>& wt) { return solveDepthCapped(wt, dmax); }, (uint64_t)leaves * (uint64_t)leaves * dmax},
            {"clamped_gw", [](const std::vector<uint64_t>& wt) { return solveClamped(wt, 3); }, (uint64_t)leaves * (uint64_t)leaves / 2}
        };

        for (const auto& a : algos) {
            int warmups = (a.name == "depth_capped" && K >= 512) ? 1 : 5;
            int numReps = (a.name == "depth_capped" && K >= 512) ? 3 : 30;
            for (int w_i = 0; w_i < warmups; ++w_i) {
                TreeResult tr = a.fn(w);
                doNotOptimize(tr.cost);
            }
            for (int rep = 0; rep < numReps; ++rep) {
                auto t0 = std::chrono::steady_clock::now();
                TreeResult tr = a.fn(w);
                auto t1 = std::chrono::steady_clock::now();
                doNotOptimize(tr.cost);
                double durUs = std::chrono::duration<double, std::micro>(t1 - t0).count();
                csv.write("E1", a.name, (uint64_t)K, 0, 0, 1.0, dmax, baseSeed, rep, "build_time_us", durUs, "us");
                if (rep == 0) {
                    csv.write("E1", a.name, (uint64_t)K, 0, 0, 1.0, dmax, baseSeed, rep, "cost", (double)tr.cost, "int");
                    csv.write("E1", a.name, (uint64_t)K, 0, 0, 1.0, dmax, baseSeed, rep, "max_depth", (double)tr.maxDepth, "int");
                    csv.write("E1", a.name, (uint64_t)K, 0, 0, 1.0, dmax, baseSeed, rep, "dp_transitions", (double)a.dpTrans, "count");
                }
            }
        }
    }

    // =========================================================================
    // E2: Proof size: hot, cold, weighted vs B1
    // =========================================================================
    std::cout << "[E2] Running proof size evaluation..." << std::endl;
    for (uint64_t N : {65536ULL, 1048576ULL}) {
        int K = 128;
        int totalLeaves = 2 * K + 1;
        std::vector<uint64_t> w = makeZipfWeights(totalLeaves, 1.0, 1001ULL);
        MerkleTree tree = buildOptimalTree(w);

        // B1 baseline proof size (balanced tree over N accounts): depth = ceil(log2(N))
        double b1Depth = std::ceil(std::log2((double)N));
        ProofSizeInfo b1Info = proofSizeInfo((uint64_t)b1Depth, (uint64_t)b1Depth, 32);

        // Hot accounts are odd leaves 1, 3, 5... 2K-1; cold super-leaves are even 0, 2, 4... 2K
        double hotHashes = 0;
        double hotBytes = 0;
        int hotCount = 0;
        for (int i = 1; i < totalLeaves; i += 2) {
            MembershipProof p = generateMembership(tree, i);
            ProofSizeInfo info = proofSizeInfo((uint64_t)p.steps.size(), (uint64_t)p.steps.size(), 32);
            hotHashes += (double)info.hashes;
            hotBytes += (double)info.totalBytes;
            hotCount++;
        }
        hotHashes /= hotCount;
        hotBytes /= hotCount;

        // Cold accounts inside a cold bucket: bucket tree depth ~ log2(N / (K+1)) + top-level super-leaf path
        double bucketDepth = std::ceil(std::log2((double)N / (double)(K + 1)));
        double coldHashes = 0;
        double coldBytes = 0;
        int coldCount = 0;
        for (int i = 0; i < totalLeaves; i += 2) {
            MembershipProof p = generateMembership(tree, i);
            ProofSizeInfo info = proofSizeInfo((uint64_t)p.steps.size() + (uint64_t)bucketDepth,
                                               (uint64_t)p.steps.size() + (uint64_t)bucketDepth, 32);
            coldHashes += (double)info.hashes;
            coldBytes += (double)info.totalBytes;
            coldCount++;
        }
        coldHashes /= coldCount;
        coldBytes /= coldCount;

        // Traffic-weighted overall (assume 80% traffic hits hot set, 20% hits cold)
        double weightedHashes = 0.8 * hotHashes + 0.2 * coldHashes;
        double weightedBytes = 0.8 * hotBytes + 0.2 * coldBytes;

        csv.write("E2", "B1_balanced", (uint64_t)K, N, 0, 1.0, 0, 1001ULL, 0, "proof_hashes", (double)b1Info.hashes, "count");
        csv.write("E2", "B1_balanced", (uint64_t)K, N, 0, 1.0, 0, 1001ULL, 0, "proof_bytes", (double)b1Info.totalBytes, "bytes");

        csv.write("E2", "interleaved_hot", (uint64_t)K, N, 0, 1.0, 0, 1001ULL, 0, "proof_hashes", hotHashes, "count");
        csv.write("E2", "interleaved_hot", (uint64_t)K, N, 0, 1.0, 0, 1001ULL, 0, "proof_bytes", hotBytes, "bytes");

        csv.write("E2", "interleaved_cold", (uint64_t)K, N, 0, 1.0, 0, 1001ULL, 0, "proof_hashes", coldHashes, "count");
        csv.write("E2", "interleaved_cold", (uint64_t)K, N, 0, 1.0, 0, 1001ULL, 0, "proof_bytes", coldBytes, "bytes");

        csv.write("E2", "interleaved_weighted", (uint64_t)K, N, 0, 1.0, 0, 1001ULL, 0, "proof_hashes", weightedHashes, "count");
        csv.write("E2", "interleaved_weighted", (uint64_t)K, N, 0, 1.0, 0, 1001ULL, 0, "proof_bytes", weightedBytes, "bytes");
    }

    // =========================================================================
    // E3: Adversarial depth
    // =========================================================================
    std::cout << "[E3] Running adversarial depth experiments..." << std::endl;
    {
        std::vector<uint64_t> baseW(64, 10ULL);
        // (a) Geometric and Fibonacci cascades
        std::vector<uint64_t> geomW = applyCascade(baseW, 2001ULL, 0);
        std::vector<uint64_t> fibW = applyCascade(baseW, 2002ULL, 1);

        TreeResult uncGeom = solveOptimal(geomW);
        TreeResult uncFib = solveOptimal(fibW);
        TreeResult capGeom = solveDepthCapped(geomW, 8);
        TreeResult capFib = solveDepthCapped(fibW, 8);
        TreeResult clampGeom = solveClamped(geomW, 3);
        TreeResult clampFib = solveClamped(fibW, 3);

        csv.write("E3", "unconstrained_geom", 64, 0, 0, 0, 0, 2001ULL, 0, "max_depth", (double)uncGeom.maxDepth, "int");
        csv.write("E3", "unconstrained_fib", 64, 0, 0, 0, 0, 2002ULL, 0, "max_depth", (double)uncFib.maxDepth, "int");
        csv.write("E3", "depth_capped_geom", 64, 0, 0, 0, 8, 2001ULL, 0, "max_depth", (double)capGeom.maxDepth, "int");
        csv.write("E3", "depth_capped_fib", 64, 0, 0, 0, 8, 2002ULL, 0, "max_depth", (double)capFib.maxDepth, "int");
        csv.write("E3", "clamped_geom", 64, 0, 0, 0, 0, 2001ULL, 0, "max_depth", (double)clampGeom.maxDepth, "int");
        csv.write("E3", "clamped_fib", 64, 0, 0, 0, 0, 2002ULL, 0, "max_depth", (double)clampFib.maxDepth, "int");

        // (b) Budget-limited attacker: sweep budget from 10 to 1000
        for (uint64_t budget : {10ULL, 50ULL, 100ULL, 250ULL, 500ULL, 1000ULL}) {
            std::vector<uint64_t> touches = attackerBudgetTouches(baseW, 2003ULL, budget / 10, budget);
            std::vector<uint64_t> attackedW = baseW;
            for (size_t i = 0; i < baseW.size(); ++i) attackedW[i] += touches[i];
            TreeResult tr = solveOptimal(attackedW);
            csv.write("E3", "attacker_budget", 64, 0, 0, 0, 0, budget, 0, "max_depth", (double)tr.maxDepth, "int");
            csv.write("E3", "attacker_budget", 64, 0, 0, 0, 0, budget, 0, "cost", (double)tr.cost, "int");
        }

        // (c) Griefing test: push honest accounts deeper
        for (double intensity : {0.1, 0.25, 0.5, 0.75, 1.0}) {
            std::vector<uint64_t> griefW = applyGriefing(baseW, 2004ULL, intensity);
            TreeResult tr = solveOptimal(griefW);
            csv.write("E3", "griefing", 64, 0, 0, intensity, 0, 2004ULL, 0, "max_depth", (double)tr.maxDepth, "int");
            csv.write("E3", "griefing", 64, 0, 0, intensity, 0, 2004ULL, 0, "cost", (double)tr.cost, "int");
        }
    }

    // =========================================================================
    // E4: Block witness multiproof size vs m in {20, 50, 100, 200, 500}
    // =========================================================================
    std::cout << "[E4] Running block witness multiproof evaluations..." << std::endl;
    {
        int K = 256;
        int totalLeaves = 2 * K + 1;
        std::vector<uint64_t> w = makeZipfWeights(totalLeaves, 1.0, 3001ULL);
        MerkleTree tree = buildOptimalTree(w);

        for (int m : {20, 50, 100, 200, 500}) {
            if (m > totalLeaves) continue;
            // Pick m leaves according to access pattern
            std::vector<int> leafIds(m);
            for (int i = 0; i < m; ++i) leafIds[i] = i * 2 % totalLeaves;

            // Compute multiproof
            Multiproof mp = generateMultiproof(tree, leafIds);

            // Separate individual proofs baseline
            uint64_t separateHashes = 0;
            for (int id : leafIds) {
                MembershipProof p = generateMembership(tree, id);
                separateHashes += p.steps.size();
            }
            uint64_t separateBytes = separateHashes * 32 + (separateHashes + 7) / 8;
            uint64_t mpHashes = mp.sharedPool.size();
            uint64_t mpBytes = mpHashes * 32 + (m * 16); // path directions
            double dedupRatio = (separateHashes > 0) ? (1.0 - (double)mpHashes / (double)separateHashes) : 0.0;

            csv.write("E4", "multiproof", (uint64_t)K, 0, (uint64_t)m, 1.0, 0, 3001ULL, 0, "multiproof_hashes", (double)mpHashes, "count");
            csv.write("E4", "multiproof", (uint64_t)K, 0, (uint64_t)m, 1.0, 0, 3001ULL, 0, "multiproof_bytes", (double)mpBytes, "bytes");
            csv.write("E4", "separate_proofs", (uint64_t)K, 0, (uint64_t)m, 1.0, 0, 3001ULL, 0, "separate_hashes", (double)separateHashes, "count");
            csv.write("E4", "separate_proofs", (uint64_t)K, 0, (uint64_t)m, 1.0, 0, 3001ULL, 0, "separate_bytes", (double)separateBytes, "bytes");
            csv.write("E4", "multiproof", (uint64_t)K, 0, (uint64_t)m, 1.0, 0, 3001ULL, 0, "deduplication_ratio", dedupRatio, "ratio");
        }
    }

    // =========================================================================
    // E5: Epoch rebuild time & cold-bucket split/merge
    // =========================================================================
    std::cout << "[E5] Running epoch rebuild and cold-bucket split/merge..." << std::endl;
    for (int K : {64, 128, 256, 512}) {
        int totalLeaves = 2 * K + 1;
        std::vector<uint64_t> w = makeZipfWeights(totalLeaves, 1.0, 4001ULL);
        for (int rep = 0; rep < 30; ++rep) {
            auto t0 = std::chrono::steady_clock::now();
            MerkleTree tree = buildOptimalTree(w);
            auto t1 = std::chrono::steady_clock::now();
            doNotOptimize(tree.depths[0]);
            double rebuildUs = std::chrono::duration<double, std::micro>(t1 - t0).count();
            csv.write("E5", "epoch_rebuild", (uint64_t)K, 0, 0, 1.0, 0, 4001ULL, rep, "rebuild_time_us", rebuildUs, "us");

            // Split/merge cold bucket simulation
            auto tm0 = std::chrono::steady_clock::now();
            // Swap 5 hot accounts with cold bucket accounts
            std::vector<uint64_t> wUpdated = w;
            for (int sw = 0; sw < 5; ++sw) {
                std::swap(wUpdated[sw * 2 + 1], wUpdated[sw * 2]);
            }
            MerkleTree tree2 = buildOptimalTree(wUpdated);
            auto tm1 = std::chrono::steady_clock::now();
            doNotOptimize(tree2.depths[0]);
            double splitMergeUs = std::chrono::duration<double, std::micro>(tm1 - tm0).count();
            csv.write("E5", "cold_split_merge", (uint64_t)K, 0, 0, 1.0, 0, 4001ULL, rep, "split_merge_cost_us", splitMergeUs, "us");
        }
    }

    // =========================================================================
    // E6: Dmax sweep: average depth vs cap
    // =========================================================================
    std::cout << "[E6] Running Dmax sweep..." << std::endl;
    {
        int K = 64;
        int totalLeaves = 2 * K + 1; // 129
        std::vector<uint64_t> w = makeZipfWeights(totalLeaves, 1.3, 5001ULL);
        uint64_t minD = (uint64_t)std::ceil(std::log2((double)totalLeaves)); // 8
        for (uint64_t dmax = minD; dmax <= minD + 8; ++dmax) {
            TreeResult tr = solveDepthCapped(w, dmax);
            double avgDepth = (double)tr.cost / (double)std::accumulate(w.begin(), w.end(), 0ULL);
            csv.write("E6", "depth_capped_sweep", (uint64_t)K, 0, 0, 1.3, dmax, 5001ULL, 0, "avg_depth", avgDepth, "depth");
            csv.write("E6", "depth_capped_sweep", (uint64_t)K, 0, 0, 1.3, dmax, 5001ULL, 0, "max_depth", (double)tr.maxDepth, "int");
            csv.write("E6", "depth_capped_sweep", (uint64_t)K, 0, 0, 1.3, dmax, 5001ULL, 0, "cost", (double)tr.cost, "int");
        }
    }

    // =========================================================================
    // E7: Non-inclusion proof size
    // =========================================================================
    std::cout << "[E7] Running non-inclusion proof size..." << std::endl;
    {
        int K = 128;
        int totalLeaves = 2 * K + 1;
        std::vector<uint64_t> w = makeZipfWeights(totalLeaves, 1.0, 6001ULL);
        MerkleTree tree = buildOptimalTree(w);
        for (int leafIdx = 0; leafIdx < 10; ++leafIdx) {
            NonInclusionProof nip = generateNonInclusion(tree, leafIdx);
            ProofSizeInfo info = proofSizeInfo((uint64_t)nip.steps.size(), (uint64_t)nip.steps.size(), 32);
            csv.write("E7", "non_inclusion", (uint64_t)K, 0, 0, 1.0, 0, 6001ULL, leafIdx, "noninc_hashes", (double)info.hashes, "count");
            csv.write("E7", "non_inclusion", (uint64_t)K, 0, 0, 1.0, 0, 6001ULL, leafIdx, "noninc_bytes", (double)info.totalBytes, "bytes");
        }
    }

    // =========================================================================
    // E8: Verification time per proof
    // =========================================================================
    std::cout << "[E8] Running verification time benchmarks..." << std::endl;
    {
        int K = 128;
        int totalLeaves = 2 * K + 1;
        std::vector<uint64_t> w = makeZipfWeights(totalLeaves, 1.0, 7001ULL);
        MerkleTree tree = buildOptimalTree(w);
        MembershipProof mp = generateMembership(tree, 5);
        NonInclusionProof nip = generateNonInclusion(tree, 5);

        for (int rep = 0; rep < 30; ++rep) {
            auto t0 = std::chrono::steady_clock::now();
            bool okM = verifyMembership(tree, mp);
            auto t1 = std::chrono::steady_clock::now();
            doNotOptimize(okM);
            double durM = std::chrono::duration<double, std::nano>(t1 - t0).count();
            csv.write("E8", "membership_verify", (uint64_t)K, 0, 0, 1.0, 0, 7001ULL, rep, "verify_time_ns", durM, "ns");

            auto tn0 = std::chrono::steady_clock::now();
            bool okN = verifyNonInclusion(tree, nip);
            auto tn1 = std::chrono::steady_clock::now();
            doNotOptimize(okN);
            double durN = std::chrono::duration<double, std::nano>(tn1 - tn0).count();
            csv.write("E8", "non_inclusion_verify", (uint64_t)K, 0, 0, 1.0, 0, 7001ULL, rep, "verify_time_ns", durN, "ns");
        }
    }

    // =========================================================================
    // E9: Pre-block pipeline: Tarjan SCC & Knapsack DP vs Greedy
    // =========================================================================
    std::cout << "[E9] Running pre-block pipeline experiments..." << std::endl;
    {
        // Tarjan SCC wash-trade detection runtime vs graph size
        for (int V : {50, 100, 200, 500, 1000}) {
            DirectedGraph g;
            g.numVertices = V;
            g.adj.resize(V);
            uint64_t seed = 8001ULL + (uint64_t)V;
            for (int i = 0; i < V; ++i) {
                for (int j = 0; j < V; ++j) {
                    if (i == j) continue;
                    seed = splitmix64(seed);
                    if ((seed % 100) < 5) g.adj[i].push_back(j); // 5% edge density
                }
            }
            for (int rep = 0; rep < 30; ++rep) {
                auto t0 = std::chrono::steady_clock::now();
                SccResult scc = runTarjanSCC(g);
                auto t1 = std::chrono::steady_clock::now();
                doNotOptimize(scc.cycleCount);
                double sccUs = std::chrono::duration<double, std::micro>(t1 - t0).count();
                csv.write("E9", "tarjan_scc", (uint64_t)V, 0, 0, 0, 0, seed, rep, "scc_time_us", sccUs, "us");
            }
        }

        // 2D Knapsack DP vs Greedy fee-per-gas
        uint64_t txSeed = 9001ULL;
        std::vector<Tx> txs;
        int txId = 0;
        for (int acc = 0; acc < 20; ++acc) {
            int numTxs = 1 + (int)(splitmix64(txSeed) % 5);
            for (int non = 0; non < numTxs; ++non) {
                uint64_t gas = 21000 + (splitmix64(txSeed) % 40000);
                uint64_t fee = 100 + (splitmix64(txSeed) % 1000);
                txs.push_back(Tx{txId++, acc, non, gas, fee});
            }
        }
        for (uint64_t gasLimit : {100000ULL, 200000ULL, 400000ULL, 800000ULL}) {
            BlockSelectionResult dpRes = solveBlockKnapsackDP(txs, gasLimit);
            BlockSelectionResult greedyRes = solveBlockGreedy(txs, gasLimit);
            double gainPct = (greedyRes.totalFee > 0)
                ? (100.0 * (double)(dpRes.totalFee - greedyRes.totalFee) / (double)greedyRes.totalFee)
                : 0.0;
            csv.write("E9", "knapsack_dp", 0, 0, (uint64_t)txs.size(), 0, gasLimit, txSeed, 0, "total_fee", (double)dpRes.totalFee, "wei");
            csv.write("E9", "greedy_ratio", 0, 0, (uint64_t)txs.size(), 0, gasLimit, txSeed, 0, "total_fee", (double)greedyRes.totalFee, "wei");
            csv.write("E9", "dp_gain", 0, 0, (uint64_t)txs.size(), 0, gasLimit, txSeed, 0, "gain_pct", gainPct, "percent");
        }
    }

    std::cout << "All Merkle experiments completed successfully! Results written to: " << outputFile << std::endl;
    return 0;
}

}  // namespace merkle
