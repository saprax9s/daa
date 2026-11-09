#include <cstdint>
#include <iostream>
#include <string>
#include <vector>

#include "automaton.h"
#include "baselines.h"
#include "bench.h"
#include "dp.h"
#include "genetic_code.h"
#include "motifs.h"
#include "prng.h"
#include "scores.h"

namespace daa {
namespace selftest {
int runAllSelfTests(int& failures);
}  // namespace selftest
}  // namespace daa

static int runToolMain(int argc, char** argv);

int main(int argc, char** argv) { return runToolMain(argc, argv); }

static int runToolMain(int argc, char** argv) {
  if (argc > 1 && std::string(argv[1]) == "--selftest") {
    int failures = 0;
    daa::selftest::runAllSelfTests(failures);
    std::cout << (failures == 0 ? "SELFTEST: ALL TESTS PASSED" : "SELFTEST: TESTS FAILED") << std::endl;
    return failures == 0 ? 0 : 1;
  }

  if (argc > 1 && std::string(argv[1]) == "--generate-instances") {
    std::string outDir = "../instances";
    if (argc > 2) outDir = argv[2];
    return daa::generateInstances(outDir);
  }

  if (argc > 1 && std::string(argv[1]) == "--bench") {
    std::string outFile = "../results/codon_results.csv";
    if (argc > 2) outFile = argv[2];
    return daa::runCodonBenchmarks(outFile);
  }

  // Default run mode: generate one instance internally with a fixed seed.
  const std::uint64_t seed = 20260911ULL;
  const int n = 6;
  daa::prng::SplitMix64 rng(seed);
  daa::dp::DpInstance inst;
  inst.n = n;
  inst.protein.resize(n);
  for (int i = 0; i < n; ++i) inst.protein[i] = (int)rng.below(20);
  inst.allowed.resize(n);
  for (int i = 0; i < n; ++i) {
    inst.allowed[i] = daa::genetic::codonsForAa(inst.protein[i]);
  }
  inst.lo = 0;
  inst.hi = 3 * n;
  inst.windowK = 0;
  inst.windowLo = 0;
  inst.windowHi = 0;

  daa::motifs::MotifSet ms;
  for (int i = 0; i < 3; ++i) {
    const int len = 2 + (int)rng.below(5);
    std::string pat;
    for (int j = 0; j < len; ++j) pat += "ACGT"[rng.below(4)];
    ms.add(pat, "random");
  }

  daa::automaton::AhoCorasick ac;
  ac.build(ms);
  const daa::scores::ScoreTable st = daa::scores::makeScoreTable(4242);
  daa::dp::DpFlags flags;
  const daa::dp::DpResult r = daa::dp::solve(ac, st, ms, inst, flags);

  std::cout << "{\"feasible\":" << (r.feasible ? "true" : "false") << ",\"score\":" << r.score
            << ",\"n\":" << n << ",\"codons\":[";
  for (std::size_t i = 0; i < r.codons.size(); ++i) {
    if (i > 0) std::cout << ",";
    std::cout << r.codons[i];
  }
  std::cout << "]}";
  std::cout << std::endl;
  return 0;
}
