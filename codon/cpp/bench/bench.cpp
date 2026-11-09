#include "bench.h"

#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <string>
#include <vector>

#include "automaton.h"
#include "baselines.h"
#include "dp.h"
#include "genetic_code.h"
#include "motifs.h"
#include "prng.h"
#include "scores.h"

namespace daa {

namespace {

volatile int64_t g_codon_sink = 0;
void doNotOptimize(int64_t v) { g_codon_sink = v; }

struct CsvWriter {
  std::ofstream out;
  explicit CsvWriter(const std::string& path) : out(path) {
    out << "language,experiment,method,n,motif_set,gc_mode,seed,rep,metric,value,unit\n";
  }

  void write(const std::string& exp, const std::string& method, int n,
             const std::string& motifSet, const std::string& gcMode,
             uint64_t seed, int rep, const std::string& metric, double val,
             const std::string& unit) {
    out << "cpp," << exp << "," << method << "," << n << "," << motifSet << ","
        << gcMode << "," << seed << "," << rep << "," << metric << ","
        << val << "," << unit << "\n";
    out.flush();
  }
};

// 15 realistic restriction enzyme recognition sites (verified against REBASE)
// EcoRI GAATTC, BamHI GGATCC, HindIII AAGCTT, XhoI CTCGAG, NdeI CATATG, NcoI CCATGG,
// XbaI TCTAGA, SpeI ACTAGT, PstI CTGCAG, KpnI GGTACC, SacI GAGCTC, NotI GCGGCCGC,
// BsaI GGTCTC, BsmBI CGTCTC, SapI GCTCTTC
motifs::MotifSet makeRealisticRestrictionMotifs() {
  motifs::MotifSet ms;
  ms.add("GAATTC", "EcoRI_verified");
  ms.add("GGATCC", "BamHI_verified");
  ms.add("AAGCTT", "HindIII_verified");
  ms.add("CTCGAG", "XhoI_verified");
  ms.add("CATATG", "NdeI_verified");
  ms.add("CCATGG", "NcoI_verified");
  ms.add("TCTAGA", "XbaI_verified");
  ms.add("ACTAGT", "SpeI_verified");
  ms.add("CTGCAG", "PstI_verified");
  ms.add("GGTACC", "KpnI_verified");
  ms.add("GAGCTC", "SacI_verified");
  ms.add("GCGGCCGC", "NotI_verified");
  ms.add("GGTCTC", "BsaI_verified");
  ms.add("CGTCTC", "BsmBI_verified");
  ms.add("GCTCTTC", "SapI_verified");
  return ms;
}

motifs::MotifSet makeRandomMotifs(int count, uint64_t seed) {
  motifs::MotifSet ms;
  prng::SplitMix64 rng(seed);
  for (int i = 0; i < count; ++i) {
    int len = 4 + (int)(rng.below(4));
    std::string pat;
    for (int j = 0; j < len; ++j) pat += "ACGT"[rng.below(4)];
    ms.add(pat, "random");
  }
  return ms;
}

// Canonical real proteins
// 1. GFP (Aequorea victoria, UniProt P42212, GenBank AAB02574.1, 238 aa)
const std::string kGfpSeq =
    "MSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTFSYGVQCF"
    "SRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLE"
    "YNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNE"
    "KRDHMVLLEFVTAAGITHGMDELYK";

// 2. Human mature insulin (UniProt P01308, B chain 30 aa + A chain 21 aa = 51 aa)
const std::string kInsulinSeq =
    "FVNQHLCGSHLVEALYLVCGERGFFYTPKTGIVEQCCTSICSLYQLENYCN";

dp::DpInstance makeInstanceFromSeq(const std::string& seq, int lo, int hi,
                                  int winK = 0, int winLo = 0, int winHi = 0,
                                  bool reduced = false) {
  dp::DpInstance inst;
  inst.n = (int)seq.size();
  inst.protein.resize(inst.n);
  inst.allowed.resize(inst.n);
  for (int i = 0; i < inst.n; ++i) {
    int aa = genetic::aaIndex(seq[i]);
    if (aa < 0) aa = 0;
    inst.protein[i] = aa;
    inst.allowed[i] = reduced ? genetic::codonsForAaReduced(aa) : genetic::codonsForAa(aa);
  }
  inst.lo = lo;
  inst.hi = hi;
  inst.windowK = winK;
  inst.windowLo = winLo;
  inst.windowHi = winHi;
  return inst;
}

dp::DpInstance makeRandomInstance(int n, uint64_t seed, int lo, int hi,
                                  int winK = 0, int winLo = 0, int winHi = 0,
                                  bool reduced = false) {
  dp::DpInstance inst;
  inst.n = n;
  inst.protein.resize(n);
  inst.allowed.resize(n);
  prng::SplitMix64 rng(seed);
  for (int i = 0; i < n; ++i) {
    int aa = (int)rng.below(20);
    inst.protein[i] = aa;
    inst.allowed[i] = reduced ? genetic::codonsForAaReduced(aa) : genetic::codonsForAa(aa);
  }
  inst.lo = lo;
  inst.hi = hi;
  inst.windowK = winK;
  inst.windowLo = winLo;
  inst.windowHi = winHi;
  return inst;
}

void writeInstanceJson(const std::string& filepath, const std::string& name,
                       const dp::DpInstance& inst, const motifs::MotifSet& ms) {
  std::ofstream out(filepath);
  out << "{\n";
  out << "  \"name\": \"" << name << "\",\n";
  out << "  \"n\": " << inst.n << ",\n";
  out << "  \"lo\": " << inst.lo << ",\n";
  out << "  \"hi\": " << inst.hi << ",\n";
  out << "  \"windowK\": " << inst.windowK << ",\n";
  out << "  \"windowLo\": " << inst.windowLo << ",\n";
  out << "  \"windowHi\": " << inst.windowHi << ",\n";
  out << "  \"protein\": [";
  for (size_t i = 0; i < inst.protein.size(); ++i) {
    if (i > 0) out << ", ";
    out << inst.protein[i];
  }
  out << "],\n";
  out << "  \"motifs\": [\n";
  for (size_t i = 0; i < ms.patterns.size(); ++i) {
    out << "    {\"pattern\": \"" << ms.patterns[i] << "\", \"source\": \""
        << ms.sources[i] << "\"}";
    if (i + 1 < ms.patterns.size()) out << ",";
    out << "\n";
  }
  out << "  ]\n";
  out << "}\n";
}

}  // namespace

int generateInstances(const std::string& instancesDir) {
  std::cout << "Generating benchmark instances in: " << instancesDir << std::endl;
  motifs::MotifSet testMs = makeRandomMotifs(3, 101ULL);
  dp::DpInstance instN5 = makeRandomInstance(5, 111ULL, 0, 15, 0, 0, 0, false);
  writeInstanceJson(instancesDir + "/small_full_n5.json", "small_full_n5", instN5, testMs);

  dp::DpInstance instN8 = makeRandomInstance(8, 222ULL, 0, 24, 0, 0, 0, true);
  writeInstanceJson(instancesDir + "/small_reduced_n8.json", "small_reduced_n8", instN8, testMs);

  dp::DpInstance instN12 = makeRandomInstance(12, 333ULL, 5, 30, 0, 0, 0, true);
  writeInstanceJson(instancesDir + "/small_reduced_n12.json", "small_reduced_n12", instN12, testMs);

  motifs::MotifSet realMs = makeRealisticRestrictionMotifs();
  dp::DpInstance gfpInst = makeInstanceFromSeq(kGfpSeq, 200, 500);
  writeInstanceJson(instancesDir + "/real_gfp.json", "real_gfp_P42212", gfpInst, realMs);

  dp::DpInstance insInst = makeInstanceFromSeq(kInsulinSeq, 50, 110);
  writeInstanceJson(instancesDir + "/real_insulin.json", "real_insulin_P01308", insInst, realMs);

  std::cout << "Instances generated successfully." << std::endl;
  return 0;
}

int runCodonBenchmarks(const std::string& outputFile) {
  CsvWriter csv(outputFile);
  std::cout << "Starting Codon Design benchmark experiments..." << std::endl;
  const uint64_t baseSeed = 20260911ULL;
  const scores::ScoreTable st = scores::makeScoreTable(baseSeed);

  // =========================================================================
  // E1: Optimality gap and validity rate: Greedy and Sample-and-Reject vs DP
  // =========================================================================
  std::cout << "[E1] Evaluating optimality gap and validity rate vs DP..." << std::endl;
  for (int seed = 1; seed <= 20; ++seed) {
    int n = 15;
    dp::DpInstance inst = makeRandomInstance(n, baseSeed + seed, 15, 30);
    motifs::MotifSet ms = makeRandomMotifs(4, 999ULL + seed);
    automaton::AhoCorasick ac;
    ac.build(ms);

    dp::DpFlags flags;
    flags.gcPruning = true;
    dp::DpResult dpRes = dp::solve(ac, st, ms, inst, flags);

    baselines::GreedyResult greedyRes = baselines::greedy(inst, st, ac, ms, true);
    baselines::SampleResult sampleRes = baselines::sampleAndReject(inst, st, ms, 5000ULL + seed, 500, 10.0);

    // Write DP ground truth
    csv.write("E1", "exact_dp", n, "random_4", "standard", seed, 0, "feasible", dpRes.feasible ? 1.0 : 0.0, "bool");
    csv.write("E1", "exact_dp", n, "random_4", "standard", seed, 0, "score", (double)dpRes.score, "score");

    // Greedy
    csv.write("E1", "greedy_m3", n, "random_4", "standard", seed, 0, "feasible", greedyRes.feasible ? 1.0 : 0.0, "bool");
    csv.write("E1", "greedy_m3", n, "random_4", "standard", seed, 0, "score", (double)greedyRes.score, "score");
    double greedyGap = (dpRes.feasible && greedyRes.feasible)
        ? (100.0 * (double)(dpRes.score - greedyRes.score) / (double)std::max((int64_t)1, std::abs(dpRes.score)))
        : 100.0;
    csv.write("E1", "greedy_m3", n, "random_4", "standard", seed, 0, "optimality_gap_pct", greedyGap, "percent");

    // Sample-and-reject
    csv.write("E1", "sample_reject_m4", n, "random_4", "standard", seed, 0, "feasible", sampleRes.feasible ? 1.0 : 0.0, "bool");
    csv.write("E1", "sample_reject_m4", n, "random_4", "standard", seed, 0, "score", (double)sampleRes.score, "score");
    double sampleGap = (dpRes.feasible && sampleRes.feasible)
        ? (100.0 * (double)(dpRes.score - sampleRes.score) / (double)std::max((int64_t)1, std::abs(dpRes.score)))
        : 100.0;
    csv.write("E1", "sample_reject_m4", n, "random_4", "standard", seed, 0, "optimality_gap_pct", sampleGap, "percent");
    csv.write("E1", "sample_reject_m4", n, "random_4", "standard", seed, 0, "wasted_candidates", (double)sampleRes.wasted, "count");
  }

  // =========================================================================
  // E2: Runtime and DP state count vs n for several motif-set sizes (|Q|)
  // =========================================================================
  std::cout << "[E2] Running scaling experiments (n vs runtime and states)..." << std::endl;
  for (int motifCount : {3, 7, 15}) {
    motifs::MotifSet ms = makeRandomMotifs(motifCount, 777ULL + motifCount);
    automaton::AhoCorasick ac;
    int qSize = ac.build(ms);
    std::string msetName = "random_" + std::to_string(motifCount);

    for (int n : {10, 20, 50, 100, 200}) {
      dp::DpInstance inst = makeRandomInstance(n, 1234ULL + n, 0, 3 * n);
      dp::DpFlags flags;
      flags.gcPruning = true;

      // Warm-up 5 runs
      for (int w = 0; w < 5; ++w) {
        dp::DpResult r = dp::solve(ac, st, ms, inst, flags);
        doNotOptimize(r.score);
      }

      // 30 repetitions
      for (int rep = 0; rep < 30; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        dp::DpResult r = dp::solve(ac, st, ms, inst, flags);
        auto t1 = std::chrono::steady_clock::now();
        doNotOptimize(r.score);
        double durUs = std::chrono::duration<double, std::micro>(t1 - t0).count();

        csv.write("E2", "exact_dp", n, msetName, "unconstrained", baseSeed, rep, "runtime_us", durUs, "us");
        if (rep == 0) {
          csv.write("E2", "exact_dp", n, msetName, "unconstrained", baseSeed, rep, "dp_states", (double)r.totalStates, "count");
          csv.write("E2", "exact_dp", n, msetName, "unconstrained", baseSeed, rep, "peak_states", (double)r.peakStates, "count");
          csv.write("E2", "exact_dp", n, msetName, "unconstrained", baseSeed, rep, "automaton_states", (double)qSize, "count");
        }
      }
    }
  }

  // =========================================================================
  // E3: Runtime and state count vs |Q| at fixed n=50
  // =========================================================================
  std::cout << "[E3] Running runtime vs |Q| at fixed n=50..." << std::endl;
  {
    int n = 50;
    dp::DpInstance inst = makeRandomInstance(n, 9876ULL, 0, 3 * n);
    for (int count : {2, 5, 10, 15, 25, 40}) {
      motifs::MotifSet ms = makeRandomMotifs(count, 5555ULL + count);
      automaton::AhoCorasick ac;
      int qSize = ac.build(ms);

      dp::DpFlags flags;
      flags.gcPruning = true;

      for (int rep = 0; rep < 30; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        dp::DpResult r = dp::solve(ac, st, ms, inst, flags);
        auto t1 = std::chrono::steady_clock::now();
        doNotOptimize(r.score);
        double durUs = std::chrono::duration<double, std::micro>(t1 - t0).count();

        csv.write("E3", "exact_dp", n, std::to_string(qSize) + "_states", "unconstrained", baseSeed, rep, "runtime_us", durUs, "us");
        if (rep == 0) {
          csv.write("E3", "exact_dp", n, std::to_string(qSize) + "_states", "unconstrained", baseSeed, rep, "dp_states", (double)r.totalStates, "count");
          csv.write("E3", "exact_dp", n, std::to_string(qSize) + "_states", "unconstrained", baseSeed, rep, "automaton_states", (double)qSize, "count");
        }
      }
    }
  }

  // =========================================================================
  // E4: Pruning ablation: D0 (none), D1 (GC bounds), D2 (dominance when lo <= 0)
  // =========================================================================
  std::cout << "[E4] Running pruning ablation..." << std::endl;
  {
    int n = 30;
    dp::DpInstance inst = makeRandomInstance(n, 3456ULL, 0, 2 * n); // lo <= 0 allows dominance
    motifs::MotifSet ms = makeRandomMotifs(5, 7890ULL);
    automaton::AhoCorasick ac;
    ac.build(ms);

    struct PruningMode {
      std::string name;
      dp::DpFlags flags;
    };
    std::vector<PruningMode> modes = {
      {"D0_no_pruning", dp::DpFlags{false, false}},
      {"D1_gc_pruning", dp::DpFlags{true, false}},
      {"D2_dominance", dp::DpFlags{true, true}}
    };

    for (const auto& pm : modes) {
      for (int rep = 0; rep < 30; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        dp::DpResult r = dp::solve(ac, st, ms, inst, pm.flags);
        auto t1 = std::chrono::steady_clock::now();
        doNotOptimize(r.score);
        double durUs = std::chrono::duration<double, std::micro>(t1 - t0).count();

        csv.write("E4", pm.name, n, "random_5", "lo_le_0", baseSeed, rep, "runtime_us", durUs, "us");
        if (rep == 0) {
          csv.write("E4", pm.name, n, "random_5", "lo_le_0", baseSeed, rep, "dp_states", (double)r.totalStates, "count");
          csv.write("E4", pm.name, n, "random_5", "lo_le_0", baseSeed, rep, "score", (double)r.score, "score");
        }
      }
    }
  }

  // =========================================================================
  // E5: Window-GC cost vs window size k
  // =========================================================================
  std::cout << "[E5] Running window-GC scaling vs window size k..." << std::endl;
  {
    int n = 20;
    motifs::MotifSet ms = makeRandomMotifs(3, 4444ULL);
    automaton::AhoCorasick ac;
    ac.build(ms);

    for (int k : {0, 2, 3, 4, 5}) {
      dp::DpInstance inst = makeRandomInstance(n, 8888ULL, 0, 3 * n, k, (k > 0 ? k : 0), (k > 0 ? 2 * k : 0));
      dp::DpFlags flags;
      flags.gcPruning = true;

      for (int rep = 0; rep < 30; ++rep) {
        auto t0 = std::chrono::steady_clock::now();
        dp::DpResult r = dp::solve(ac, st, ms, inst, flags);
        auto t1 = std::chrono::steady_clock::now();
        doNotOptimize(r.score);
        double durUs = std::chrono::duration<double, std::micro>(t1 - t0).count();

        csv.write("E5", "window_gc_k" + std::to_string(k), n, "random_3", "window", baseSeed, rep, "runtime_us", durUs, "us");
        if (rep == 0) {
          csv.write("E5", "window_gc_k" + std::to_string(k), n, "random_3", "window", baseSeed, rep, "dp_states", (double)r.totalStates, "count");
          csv.write("E5", "window_gc_k" + std::to_string(k), n, "random_3", "window", baseSeed, rep, "feasible", r.feasible ? 1.0 : 0.0, "bool");
        }
      }
    }
  }

  // =========================================================================
  // E6: Realistic restriction-site set on random and real proteins (GFP / Insulin)
  // =========================================================================
  std::cout << "[E6] Running realistic restriction sites on GFP and Insulin..." << std::endl;
  {
    motifs::MotifSet realMs = makeRealisticRestrictionMotifs();
    automaton::AhoCorasick ac;
    int qSize = ac.build(realMs);

    struct TargetProtein {
      std::string name;
      std::string seq;
      int tightLo;
      int tightHi;
      int looseLo;
      int looseHi;
    };

    std::vector<TargetProtein> targets = {
      {"insulin_51aa", kInsulinSeq, 60, 90, 30, 120},
      {"gfp_238aa", kGfpSeq, 280, 350, 150, 500}
    };

    for (const auto& tp : targets) {
      for (const std::string gcSetting : {"tight", "loose"}) {
        int lo = (gcSetting == "tight") ? tp.tightLo : tp.looseLo;
        int hi = (gcSetting == "tight") ? tp.tightHi : tp.looseHi;
        dp::DpInstance inst = makeInstanceFromSeq(tp.seq, lo, hi);

        dp::DpFlags flags;
        flags.gcPruning = true;

        for (int rep = 0; rep < 30; ++rep) {
          auto t0 = std::chrono::steady_clock::now();
          dp::DpResult r = dp::solve(ac, st, realMs, inst, flags);
          auto t1 = std::chrono::steady_clock::now();
          doNotOptimize(r.score);
          double durUs = std::chrono::duration<double, std::micro>(t1 - t0).count();

          csv.write("E6", tp.name, inst.n, "realistic_15_enzymes", gcSetting, baseSeed, rep, "runtime_us", durUs, "us");
          if (rep == 0) {
            csv.write("E6", tp.name, inst.n, "realistic_15_enzymes", gcSetting, baseSeed, rep, "feasible", r.feasible ? 1.0 : 0.0, "bool");
            csv.write("E6", tp.name, inst.n, "realistic_15_enzymes", gcSetting, baseSeed, rep, "score", (double)r.score, "score");
            csv.write("E6", tp.name, inst.n, "realistic_15_enzymes", gcSetting, baseSeed, rep, "dp_states", (double)r.totalStates, "count");
            csv.write("E6", tp.name, inst.n, "realistic_15_enzymes", gcSetting, baseSeed, rep, "automaton_states", (double)qSize, "count");
          }
        }
      }
    }
  }

  // =========================================================================
  // E7: Wasted candidates of sample-and-reject vs constraint tightness
  // =========================================================================
  std::cout << "[E7] Running wasted candidates vs constraint tightness..." << std::endl;
  {
    int n = 20;
    motifs::MotifSet ms = makeRealisticRestrictionMotifs();
    int mid = 30;
    for (int margin : {3, 6, 12, 18, 25}) {
      dp::DpInstance inst = makeRandomInstance(n, 7777ULL, mid - margin, mid + margin);
      baselines::SampleResult sr = baselines::sampleAndReject(inst, st, ms, 9999ULL + margin, 1000, 10.0);
      double wastedPct = (sr.trials > 0) ? (100.0 * (double)sr.wasted / (double)sr.trials) : 0.0;
      csv.write("E7", "sample_and_reject", n, "realistic_15", "margin_" + std::to_string(margin),
                9999ULL + margin, 0, "wasted_candidates", (double)sr.wasted, "count");
      csv.write("E7", "sample_and_reject", n, "realistic_15", "margin_" + std::to_string(margin),
                9999ULL + margin, 0, "total_trials", (double)sr.trials, "count");
      csv.write("E7", "sample_and_reject", n, "realistic_15", "margin_" + std::to_string(margin),
                9999ULL + margin, 0, "wasted_pct", wastedPct, "percent");
    }
  }

  // =========================================================================
  // E8: Peak memory
  // =========================================================================
  std::cout << "[E8] Estimating peak memory across scaling sizes..." << std::endl;
  {
    motifs::MotifSet ms = makeRealisticRestrictionMotifs();
    automaton::AhoCorasick ac;
    ac.build(ms);
    for (int n : {20, 50, 100, 200}) {
      dp::DpInstance inst = makeRandomInstance(n, 1111ULL + n, 0, 3 * n);
      dp::DpFlags flags;
      flags.gcPruning = true;
      dp::DpResult r = dp::solve(ac, st, ms, inst, flags);
      // Entry struct is ~48 bytes on 64-bit systems
      double peakKb = (double)(r.peakStates * 48ULL) / 1024.0;
      csv.write("E8", "exact_dp", n, "realistic_15", "unconstrained", baseSeed, 0, "peak_memory_kb", peakKb, "kb");
      csv.write("E8", "exact_dp", n, "realistic_15", "unconstrained", baseSeed, 0, "peak_states", (double)r.peakStates, "count");
    }
  }

  std::cout << "All Codon benchmarks finished! Results saved to: " << outputFile << std::endl;
  return 0;
}

}  // namespace daa

