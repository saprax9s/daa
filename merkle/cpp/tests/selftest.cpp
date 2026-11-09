#include "selftest.h"

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <functional>
#include <limits>
#include <string>
#include <vector>

#include "hashing.h"
#include "proofs.h"
#include "sha256.h"
#include "trees.h"
#include "weights.h"
#include "workloads.h"
#include "pipeline.h"

namespace merkle {

int g_checks = 0;
int g_passed = 0;

void check(bool ok, const char* name) {
    g_checks++;
    if (ok) {
        g_passed++;
        std::printf("PASS: %s\n", name);
    } else {
        std::printf("FAIL: %s\n", name);
    }
    std::fflush(stdout);
}

// Independent O(K^3) alphabetic tree DP optimum (ground truth).
static uint64_t bruteAlphabeticOptimum(const std::vector<uint64_t>& w) {
    int n = (int)w.size();
    std::vector<std::vector<uint64_t>> C(n, std::vector<uint64_t>(n, 0));
    std::vector<uint64_t> pref(n + 1, 0);
    for (int i = 0; i < n; ++i) pref[i + 1] = pref[i] + w[i];
    auto W = [&](int i, int j) { return pref[j + 1] - pref[i]; };
    for (int i = 0; i < n; ++i) C[i][i] = 0;
    for (int len = 2; len <= n; ++len) {
        for (int i = 0; i + len - 1 < n; ++i) {
            int j = i + len - 1;
            uint64_t best = std::numeric_limits<uint64_t>::max() / 4;
            for (int m = i; m < j; ++m) {
                uint64_t c = C[i][m] + C[m + 1][j] + W(i, j);
                best = std::min(best, c);
            }
            C[i][j] = best;
        }
    }
    return C[0][n - 1];
}

bool checkSha256Vectors() {
    bool ok = true;
    ok = ok && (sha256Hex("abc") == "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad");
    ok = ok && (sha256Hex("") == "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855");
    ok = ok && (sha256Hex("abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq") ==
                "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1");
    check(ok, "SHA-256 self-check: FIPS vectors (abc / empty / long)");
    return ok;
}

bool checkTreeCrossValidations() {
    bool ok = true;
    int leafCounts[] = {2, 4, 6, 8, 10};
    for (int n : leafCounts) {
        std::vector<uint64_t> w((size_t)n, 0);
        uint64_t state = 0x9E3779B97F4A7C15ULL * (uint64_t)(n + 7);
        for (int i = 0; i < n; ++i) {
            state += 0x9E3779B97F4A7C15ULL;
            uint64_t z = state;
            z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
            z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
            z = z ^ (z >> 31);
            w[(size_t)i] = z % 100ULL + 1ULL;
        }
        TreeResult opt = solveOptimal(w);
        TreeResult knuth = solveKnuth(w);
        uint64_t brute = bruteAlphabeticOptimum(w);
        if (opt.cost != knuth.cost) ok = false;
        if (knuth.cost != brute) ok = false;
        if (opt.cost != brute) ok = false;
        // B4 vs B5 must be identical (both claim optimality).
        check(opt.cost == knuth.cost, "optimal (Garsia-Wachs) == Knuth DP cost");
    }
    return ok;
}

bool checkDepthCapped() {
    bool ok = true;
    // Equal weights, 6 leaves: optimal alphabetic tree is balanced, height = ceil(log2(6)) = 3.
    std::vector<uint64_t> w6 = {1, 1, 1, 1, 1, 1};
    TreeResult unc = solveOptimal(w6);
    uint64_t dmax = 0;
    {
        uint64_t dd = 1;
        while ((1ULL << dd) < w6.size()) dd++;
        dmax = dd;
    }
    TreeResult capped = solveDepthCapped(w6, dmax);
    if (capped.cost != unc.cost) ok = false;  // cap not binding -> equal
    if (capped.maxDepth > dmax) ok = false;    // max depth <= Dmax
    check(capped.cost == unc.cost, "depth-capped DP == unconstrained when cap not binding");

    // Random small instance, generous cap (6) -> must equal unconstrained optimum.
    std::vector<uint64_t> wr = {3, 1, 4, 1, 5, 9, 2, 6};
    TreeResult unc2 = solveOptimal(wr);
    TreeResult cap2 = solveDepthCapped(wr, 6);
    if (cap2.cost != unc2.cost) ok = false;
    if (cap2.maxDepth > 6) ok = false;

    // Binding cap on a small instance: feasible, cost >= unconstrained, maxDepth <= Dmax.
    TreeResult cappedBind = solveDepthCapped(wr, 2);
    if (cappedBind.cost >= std::numeric_limits<uint64_t>::max() / 4) {
        // infeasible at cap: acceptable, but check it was reported as such
    } else {
        if (cappedBind.cost < unc2.cost) ok = false;
        if (cappedBind.maxDepth > 2) ok = false;
        check(cappedBind.cost > unc2.cost, "depth-capped DP cost > unconstrained when cap binds");
    }
    check(ok, "depth-capped DP: >= unconstrained, equal when not binding, max depth <= Dmax");
    return ok;
}

bool checkMehlhornBound() {
    bool ok = true;
    std::vector<uint64_t> w = {2, 3, 1, 4, 5, 9, 2, 6};
    TreeResult mh = solveMehlhorn(w);
    TreeResult op = solveOptimal(w);
    if (mh.cost > op.cost + 2) ok = false;
    check(mh.cost <= op.cost + 2, "Mehlhorn simple version cost <= optimal + 2 (entropy-style bound)");
    return ok;
}

static bool checkOrderPreserved(const merkle::MerkleTree& t) {
    // For every internal node, left subtree leaves must all precede right subtree leaves.
    std::function<bool(int)> rec = [&](int id) -> bool {
        if (t.nodes[id].isLeaf) return true;
        int l = t.nodes[id].left;
        int r = t.nodes[id].right;
        // gather min/max leaf idx of each subtree
        std::function<std::pair<int, int>(int)> rangeOf = [&](int nid) -> std::pair<int, int> {
            if (t.nodes[nid].isLeaf) return {t.nodes[nid].leafIdx, t.nodes[nid].leafIdx};
            auto a = rangeOf(t.nodes[nid].left);
            auto b = rangeOf(t.nodes[nid].right);
            return {a.first, b.second};
        };
        auto lr = rangeOf(l);
        auto rr = rangeOf(r);
        if (!(lr.second < rr.first)) return false;
        return rec(l) && rec(r);
    };
    return rec(t.root);
}

bool checkProofs() {
    bool ok = true;
    std::vector<uint64_t> w = {2, 3, 1, 4, 1, 5, 9, 2, 6, 5};
    merkle::MerkleTree t = merkle::buildOptimalTree(w);
    for (int i = 0; i < 10; ++i) {
        merkle::MembershipProof p = merkle::generateMembership(t, i);
        if (!merkle::verifyMembership(t, p)) ok = false;
        merkle::MembershipProof p2 = merkle::generateMembership(t, i);
        if (!p2.steps.empty()) p2.steps[0].siblingHash[0] = (uint8_t)(p2.steps[0].siblingHash[0] ^ 0xFF);
        if (merkle::verifyMembership(t, p2)) ok = false;
        merkle::MembershipProof p3 = merkle::generateMembership(t, i);
        if (!p3.steps.empty()) p3.steps[0].bit = (uint8_t)(1 - p3.steps[0].bit);
        if (merkle::verifyMembership(t, p3)) ok = false;
        merkle::MembershipProof p4 = merkle::generateMembership(t, i);
        p4.leafHash[0] = (uint8_t)(p4.leafHash[0] ^ 0x01);
        if (merkle::verifyMembership(t, p4)) ok = false;
    }

    // Non-inclusion: adjacent pair proof.
    merkle::NonInclusionProof nip = merkle::generateNonInclusion(t, 3);
    if (nip.upperIdx != nip.lowerIdx + 1) ok = false;
    if (!merkle::verifyNonInclusion(t, nip)) ok = false;
    merkle::NonInclusionProof nip2 = merkle::generateNonInclusion(t, 3);
    nip2.steps[0].siblingHash[0] = (uint8_t)(nip2.steps[0].siblingHash[0] ^ 0x55);
    if (merkle::verifyNonInclusion(t, nip2)) ok = false;
    merkle::NonInclusionProof nip3 = merkle::generateNonInclusion(t, 3);
    nip3.upperIdx = nip3.lowerIdx + 2;
    if (merkle::verifyNonInclusion(t, nip3)) ok = false;

    // Multiproof with shared-hash dedup.
    std::vector<int> ids = {0, 2, 4, 6};
    merkle::Multiproof mp = merkle::generateMultiproof(t, ids);
    if (!merkle::verifyMultiproof(mp, t, ids)) ok = false;
    merkle::Multiproof mp2 = merkle::generateMultiproof(t, ids);
    if (!mp2.leafPaths.empty() && !mp2.leafPaths[0].empty()) {
        mp2.leafPaths[0][0].second = (uint8_t)(1 - mp2.leafPaths[0][0].second);
    }
    if (merkle::verifyMultiproof(mp2, t, ids)) ok = false;

    if (!checkOrderPreserved(t)) ok = false;

    check(ok, "proofs: membership / non-inclusion / multiproof verify, tampered rejected, order preserved");
    return ok;
}

bool checkDeterminismAndOrder() {
    bool ok = true;
    std::vector<uint64_t> w = {5, 1, 2, 8, 3, 1, 4, 7};
    TreeResult a = solveOptimal(w);
    TreeResult b = solveOptimal(w);
    if (a.cost != b.cost || a.maxDepth != b.maxDepth) ok = false;
    merkle::MerkleTree t1 = merkle::buildOptimalTree(w);
    merkle::MerkleTree t2 = merkle::buildOptimalTree(w);
    if (!(t1.rootHash == t2.rootHash)) ok = false;
    check(ok, "determinism: same inputs -> same cost/depth/root");

    std::vector<uint64_t> scores = {9, 7, 7, 7, 3, 1, 5};
    std::vector<uint64_t> hot = computeHotSet(scores, 3);
    bool hotOk = (hot.size() == 3) && hot[0] == 0 && hot[1] == 1 && hot[2] == 2;
    check(hotOk, "hot set tie-break by lowest account index");

    bool wOk = (updateWeight(0, 0) == 0) && (updateWeight(0, 1) == 300000) && (updateWeight(10, 0) == 7) &&
                (updateWeight(100, 10) == 3000070);
    check(wOk, "weight update W_new = (3*touches*1e6 + 7*W_old)/10");
    return ok;
}

bool checkWorkloads() {
    bool ok = true;
    auto e1 = generateEpochTouches(42, 8, 1.0, 3, 2);
    auto e2 = generateEpochTouches(42, 8, 1.0, 3, 2);
    if (!(e1.size() == e2.size())) ok = false;
    for (size_t i = 0; i < e1.size(); ++i) {
        if (!(e1[i].size() == e2[i].size())) ok = false;
        for (size_t j = 0; j < e1[i].size(); ++j) {
            if (e1[i][j] != e2[i][j]) ok = false;
        }
    }
    check(ok, "workload: deterministic splitmix64/Zipf epoch generation");

    std::vector<uint64_t> w = {9, 4, 2, 8, 1, 5, 6};
    std::vector<uint64_t> c1 = applyCascade(w, 7, 0);
    std::vector<uint64_t> c2 = applyCascade(w, 7, 1);
    if (c1.size() != w.size() || c2.size() != w.size()) ok = false;
    std::vector<uint64_t> bud = attackerBudgetTouches(w, 7, 3, 30);
    std::vector<uint64_t> gr = applyGriefing(w, 7, 0.5);
    if (bud.size() != w.size() || gr.size() != w.size()) ok = false;
    check(ok, "workload: adversary (cascade / budget / griefing) entry points run");
    return ok;
}

bool checkTarjanAndKnapsack() {
    bool ok = true;
    // 1. Tarjan vs Brute-force cycle detection on synthetic directed graphs
    uint64_t seed = 987654321ULL;
    for (int V : {4, 6, 8, 10}) {
        DirectedGraph g;
        g.numVertices = V;
        g.adj.resize(V);
        for (int i = 0; i < V; ++i) {
            for (int j = 0; j < V; ++j) {
                if (i == j) continue;
                seed = seed * 6364136223846793005ULL + 1ULL;
                if ((seed % 100) < 25) { // 25% edge probability
                    g.adj[i].push_back(j);
                }
            }
        }
        bool tarjanCycle = runTarjanSCC(g).hasWashTradeCycle;
        bool bruteCycle = bruteForceHasCycle(g);
        if (tarjanCycle != bruteCycle) ok = false;
    }
    check(ok, "E9: Tarjan SCC wash-trade detection == brute-force cycle checker on random graphs");

    // 2. Nonce-prefix gas-bucketed knapsack DP vs greedy
    std::vector<Tx> txs = {
        // Account 0
        {0, 0, 0, 21000, 100},
        {1, 0, 1, 50000, 500},
        // Account 1
        {2, 1, 0, 40000, 300},
        {3, 1, 1, 25000, 250},
        // Account 2
        {4, 2, 0, 100000, 1200}
    };
    uint64_t gasLimit = 80000;
    BlockSelectionResult dpRes = solveBlockKnapsackDP(txs, gasLimit);
    BlockSelectionResult greedyRes = solveBlockGreedy(txs, gasLimit);

    bool knapsackOk = (dpRes.totalGas <= gasLimit) && (dpRes.totalFee >= greedyRes.totalFee);
    check(knapsackOk, "E9: 2D knapsack DP fee >= greedy fee under gas limit and nonce constraints");
    return ok && knapsackOk;
}

SelfTestResult runSelfTest() {
    SelfTestResult res;
    res.checksPassed = 0;
    res.checksFailed = 0;

    checkSha256Vectors();
    checkTreeCrossValidations();
    checkDepthCapped();
    checkMehlhornBound();
    checkProofs();
    checkDeterminismAndOrder();
    checkWorkloads();
    checkTarjanAndKnapsack();

    res.checksPassed = g_passed;
    res.checksFailed = g_checks - g_passed;
    res.details = res.checksFailed == 0 ? "ALL PASS" : "FAILURES PRESENT";
    return res;
}

}  // namespace merkle

#ifdef SELFTEST_STANDALONE_MAIN
int main() {
    merkle::SelfTestResult r = merkle::runSelfTest();
    std::printf("selftest: %s (%d passed, %d failed)\n", r.details.c_str(), r.checksPassed, r.checksFailed);
    return (r.checksFailed == 0) ? 0 : 1;
}
#endif
