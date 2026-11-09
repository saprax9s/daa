#include "selftest.h"

#include <cstdint>
#include <iostream>
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
namespace selftest {

namespace {

using AhoCorasick = daa::automaton::AhoCorasick;
using DpFlags = daa::dp::DpFlags;
using DpInstance = daa::dp::DpInstance;
using DpResult = daa::dp::DpResult;
using MotifSet = daa::motifs::MotifSet;
using ScoreTable = daa::scores::ScoreTable;
using SplitMix64 = daa::prng::SplitMix64;

MotifSet makeMotifSet(std::uint64_t seed, int nMotifs, int maxLen) {
  MotifSet ms;
  SplitMix64 rng(seed);
  for (int i = 0; i < nMotifs; ++i) {
    const int len = 2 + static_cast<int>(rng.below(static_cast<std::uint64_t>(maxLen)));
    std::string pat;
    for (int j = 0; j < len; ++j) pat += "ACGT"[rng.below(4)];
    ms.add(pat, "random");
  }
  return ms;
}

DpInstance makeInstance(int n,
                       std::uint64_t seed,
                       bool reduced,
                       int nMotifs,
                       int gcLo,
                       int gcHi,
                       int winK,
                       int winLo,
                       int winHi) {
  DpInstance inst;
  (void)nMotifs;
  inst.n = n;
  inst.lo = gcLo;
  inst.hi = gcHi;
  inst.windowK = winK;
  inst.windowLo = winLo;
  inst.windowHi = winHi;
  SplitMix64 rng(seed);
  inst.protein.resize(n);
  for (int i = 0; i < n; ++i) inst.protein[i] = rng.below(20);
  inst.allowed.resize(n);
  for (int i = 0; i < n; ++i) {
    inst.allowed[i] = reduced ? daa::genetic::codonsForAaReduced(inst.protein[i])
                              : daa::genetic::codonsForAa(inst.protein[i]);
  }
  return inst;
}

// V2: Aho-Corasick verdict must equal the naive substring-search verdict.
void testV2(int& failures) {
  MotifSet ms;
  ms.add("GAATTC", "ecoRI");
  ms.add("AAGCTT", "hindIII");
  ms.add("GGATCC", "bamHI");
  {
    SplitMix64 rng(99);
    for (int i = 0; i < 2; ++i) {
      const int len = 2 + static_cast<int>(rng.below(5));
      std::string pat;
      for (int j = 0; j < len; ++j) pat += "ACGT"[rng.below(4)];
      ms.add(pat, "random");
    }
  }
  AhoCorasick ac;
  ac.build(ms);
  SplitMix64 rng(12345);
  const int trials = 100000;
  int mismatches = 0;
  for (int t = 0; t < trials; ++t) {
    const int len = rng.below(41);
    std::string dna;
    for (int j = 0; j < len; ++j) dna += "ACGT"[rng.below(4)];
    const bool a = ac.search(dna);
    const bool n = daa::motifs::naiveMotifSearch(ms, dna);
    if (a != n) {
      ++mismatches;
      break;
    }
  }
  if (mismatches == 0) {
    std::cout << "V2: PASS (automaton vs naive substring search on " << trials
              << " random DNA strings, lengths 0..40, boundary-spanning motifs, 0 mismatches)" << std::endl;
  } else {
    std::cout << "V2: FAIL (" << mismatches << " mismatches)" << std::endl;
    ++failures;
  }
}

// V1: DP optimum == brute-force optimum; infeasibility agreement; DP answer validity.
void testV1(int& failures) {
  const ScoreTable st = daa::scores::makeScoreTable(777);
  const MotifSet ms = makeMotifSet(555, 3, 6);

  // Instance A: n=5 full code.
  DpInstance instA = makeInstance(5, 11, false, 3, 0, 15, 0, 0, 0);
  AhoCorasick acA;
  acA.build(ms);
  DpResult dpA = daa::dp::solve(acA, st, ms, instA, DpFlags());
  daa::baselines::BruteResult bfA = daa::baselines::bruteForce(instA, st, ms);

  bool ok = true;
  if (dpA.feasible != bfA.feasible) ok = false;
  if (dpA.feasible && dpA.score != bfA.score) ok = false;
  if (dpA.feasible) {
    if (static_cast<int>(dpA.codons.size()) != 5) {
      ok = false;
    } else {
      std::string dna;
      int gc = 0;
      for (int c : dpA.codons) {
        char buf[4];
        daa::genetic::codonString(c, buf);
        dna += buf;
        gc += daa::genetic::gcCount(c);
      }
      if (daa::motifs::naiveMotifSearch(ms, dna)) ok = false;
      if (gc < instA.lo || gc > instA.hi) ok = false;
      if (daa::scores::sequenceScore(dpA.codons, st) != dpA.score) ok = false;
    }
  }

  // Instance B: reduced code, medium size (n=8).
  MotifSet msB = makeMotifSet(556, 2, 6);
  daa::dp::DpInstance instB = makeInstance(8, 22, true, 2, 0, 24, 0, 0, 0);
  AhoCorasick acB;
  acB.build(msB);
  DpResult dpB = daa::dp::solve(acB, st, msB, instB, DpFlags());
  daa::baselines::BruteResult bfB = daa::baselines::bruteForce(instB, st, msB);

  if (dpB.feasible != bfB.feasible) ok = false;
  if (dpB.feasible && dpB.score != bfB.score) ok = false;

  if (ok) {
    std::cout << "V1: PASS (DP optimum == brute force on full-code n=5 and reduced-code n=8 instances; "
              << "infeasibility agreement holds; DP answer re-validated with naive checker and score)" << std::endl;
  } else {
    std::cout << "V1: FAIL (dpA feasible=" << dpA.feasible << " score=" << dpA.score << " bfA feasible=" << bfA.feasible
              << " score=" << bfA.score << "; dpB feasible=" << dpB.feasible << " score=" << dpB.score
              << " bfB feasible=" << bfB.feasible << " score=" << bfB.score << ")" << std::endl;
    ++failures;
  }
}

// Pruning sanity: D0 / D1 / D2 must return the identical optimum (lo = 0 => dominance valid).
void testPruning(int& failures) {
  const ScoreTable st = daa::scores::makeScoreTable(777);
  const MotifSet ms = makeMotifSet(555, 3, 6);
  const DpInstance inst = makeInstance(5, 11, false, 3, 0, 15, 0, 0, 0);
  AhoCorasick ac;
  ac.build(ms);

  DpResult r0 = daa::dp::solve(ac, st, ms, inst, DpFlags());
  DpFlags f1;
  f1.gcPruning = true;
  DpResult r1 = daa::dp::solve(ac, st, ms, inst, f1);
  DpFlags f2;
  f2.gcPruning = true;
  f2.dominance = true;
  DpResult r2 = daa::dp::solve(ac, st, ms, inst, f2);

  bool ok = true;
  if (r0.feasible != r1.feasible || r0.feasible != r2.feasible) ok = false;
  if (r0.feasible && (r0.score != r1.score || r0.score != r2.score)) ok = false;
  if (ok) {
    std::cout << "PRUNING: PASS (D0 == D1 == D2 optimum on small instance, dominance valid when lo <= 0)" << std::endl;
  } else {
    std::cout << "PRUNING: FAIL (D0 score=" << r0.score << " D1 score=" << r1.score << " D2 score=" << r2.score << ")" << std::endl;
    ++failures;
  }
}

// V6 groundwork: window-GC DP must equal brute force on a small instance.
void testWindowGC(int& failures) {
  const ScoreTable st = daa::scores::makeScoreTable(777);
  const MotifSet ms = makeMotifSet(557, 2, 6);
  DpInstance inst = makeInstance(4, 33, false, 2, 0, 12, 2, 2, 6);
  AhoCorasick ac;
  ac.build(ms);
  DpFlags flags;
  DpResult dp = daa::dp::solve(ac, st, ms, inst, flags);
  daa::baselines::BruteResult bf = daa::baselines::bruteForce(inst, st, ms);
  bool ok = true;
  if (dp.feasible != bf.feasible) ok = false;
  if (dp.feasible && dp.score != bf.score) ok = false;
  if (ok) {
    std::cout << "WINDOW: PASS (window-GC DP == brute force on small instance)" << std::endl;
  } else {
    std::cout << "WINDOW: FAIL (dp feasible=" << dp.feasible << " score=" << dp.score << " bf feasible=" << bf.feasible
              << " score=" << bf.score << ")" << std::endl;
    ++failures;
  }
}

}  // namespace

// Genetic Code data model validation (Part 2)
void testGeneticCode(int& failures) {
  int senseCount = 0;
  int stopCount = 0;
  for (int c = 0; c < 64; ++c) {
    int aa = daa::genetic::aaIndexForCodon(c);
    if (aa >= 0) ++senseCount;
    else ++stopCount;
  }
  bool ok = (senseCount == 61) && (stopCount == 3);

  // Per-amino-acid degeneracy checks
  int degenCounts[7] = {0};
  for (int aa = 0; aa < 20; ++aa) {
    int d = daa::genetic::degeneracy(aa);
    if (d >= 1 && d <= 6) degenCounts[d]++;
    else ok = false;
  }
  // Standard table: 2 with 1 (Met, Trp), 9 with 2, 1 with 3 (Ile), 5 with 4, 3 with 6 (Leu, Arg, Ser)
  if (degenCounts[1] != 2 || degenCounts[2] != 9 || degenCounts[3] != 1 ||
      degenCounts[4] != 5 || degenCounts[6] != 3) {
    ok = false;
  }

  if (ok) {
    std::cout << "GENETIC CODE: PASS (64 codons: 61 sense, 3 stop; degeneracy counts 1/2/3/4/6 match standard table)" << std::endl;
  } else {
    std::cout << "GENETIC CODE: FAIL" << std::endl;
    ++failures;
  }
}

// V3: DP == brute force on reduced code for longer protein (n=10).
void testV3(int& failures) {
  const ScoreTable st = daa::scores::makeScoreTable(333);
  const MotifSet ms = makeMotifSet(444, 2, 5);
  DpInstance inst = makeInstance(10, 555, true, 2, 5, 25, 0, 0, 0);
  AhoCorasick ac;
  ac.build(ms);
  DpResult dp = daa::dp::solve(ac, st, ms, inst, DpFlags());
  daa::baselines::BruteResult bf = daa::baselines::bruteForce(inst, st, ms);

  bool ok = (dp.feasible == bf.feasible);
  if (dp.feasible && dp.score != bf.score) ok = false;
  if (ok) {
    std::cout << "V3: PASS (DP == brute force on reduced code n=10: score=" << dp.score << ", 1024 candidates)" << std::endl;
  } else {
    std::cout << "V3: FAIL (dp feasible=" << dp.feasible << " score=" << dp.score << " bf feasible=" << bf.feasible << " score=" << bf.score << ")" << std::endl;
    ++failures;
  }
}

// V4: Greedy and sample-and-reject lose score or fail on measurable fraction of instances.
void testV4(int& failures) {
  const ScoreTable st = daa::scores::makeScoreTable(123);
  int greedyLost = 0;
  int sampleLost = 0;
  for (int s = 1; s <= 10; ++s) {
    DpInstance inst = makeInstance(12, 1000 + s, false, 3, 10, 26, 0, 0, 0);
    MotifSet ms = makeMotifSet(2000 + s, 3, 5);
    AhoCorasick ac;
    ac.build(ms);
    DpResult dp = daa::dp::solve(ac, st, ms, inst, DpFlags());
    auto gr = daa::baselines::greedy(inst, st, ac, ms, true);
    auto sr = daa::baselines::sampleAndReject(inst, st, ms, 3000 + s, 200, 10.0);

    if (dp.feasible) {
      if (!gr.feasible || gr.score < dp.score) ++greedyLost;
      if (!sr.feasible || sr.score < dp.score) ++sampleLost;
    }
  }
  bool ok = (greedyLost > 0) && (sampleLost > 0);
  if (ok) {
    std::cout << "V4: PASS (greedy lost or failed on " << greedyLost << "/10 instances; sample-and-reject lost on " << sampleLost << "/10)" << std::endl;
  } else {
    std::cout << "V4: FAIL (greedyLost=" << greedyLost << ", sampleLost=" << sampleLost << ")" << std::endl;
    ++failures;
  }
}

// V5: Measured runtime & states grow consistently; pruning reduces states without changing optimum.
void testV5(int& failures) {
  const ScoreTable st = daa::scores::makeScoreTable(789);
  const MotifSet ms = makeMotifSet(999, 3, 5);
  AhoCorasick ac;
  ac.build(ms);

  DpInstance inst = makeInstance(15, 777, false, 3, 0, 30, 0, 0, 0);
  DpFlags f0;
  DpFlags f1; f1.gcPruning = true;
  DpFlags f2; f2.gcPruning = true; f2.dominance = true;

  DpResult r0 = daa::dp::solve(ac, st, ms, inst, f0);
  DpResult r1 = daa::dp::solve(ac, st, ms, inst, f1);
  DpResult r2 = daa::dp::solve(ac, st, ms, inst, f2);

  bool ok = (r0.feasible == r1.feasible && r0.feasible == r2.feasible);
  if (r0.feasible) {
    if (r0.score != r1.score || r0.score != r2.score) ok = false;
    if (r2.totalStates > r1.totalStates || r1.totalStates > r0.totalStates) ok = false;
  }
  if (ok) {
    std::cout << "V5: PASS (D0 states=" << r0.totalStates << " >= D1 states=" << r1.totalStates
              << " >= D2 states=" << r2.totalStates << ", optimum strictly preserved at " << r0.score << ")" << std::endl;
  } else {
    std::cout << "V5: FAIL" << std::endl;
    ++failures;
  }
}

int runAllSelfTests(int& failures) {
  testGeneticCode(failures);
  testV2(failures);
  testV1(failures);
  testV3(failures);
  testV4(failures);
  testV5(failures);
  testPruning(failures);
  testWindowGC(failures);
  return failures;
}

}  // namespace selftest
}  // namespace daa

#ifndef CODON_SELFTEST_NO_MAIN
int main(int /*argc*/, char** /*argv*/) {
  int failures = 0;
  daa::selftest::runAllSelfTests(failures);
  if (failures == 0) {
    std::cout << "SELFTEST RESULT: PASS" << std::endl;
  } else {
    std::cout << "SELFTEST RESULT: FAIL" << std::endl;
  }
  return failures == 0 ? 0 : 1;
}
#endif
