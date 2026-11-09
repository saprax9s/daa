#pragma once

#include <cstdint>
#include <vector>

#include "dp.h"
#include "genetic_code.h"
#include "motifs.h"
#include "prng.h"
#include "scores.h"

namespace daa {
namespace baselines {

struct BruteResult {
  bool feasible;
  int64_t score;
  std::vector<int> codons;
  long long evaluated;
};

// M1: brute force / backtracking ground truth. Uses the independent naive
// substring validity check, the shared score function and the GC/window rules.
BruteResult bruteForce(const daa::dp::DpInstance& inst,
                       const daa::scores::ScoreTable& st,
                       const daa::motifs::MotifSet& ms);

struct GreedyResult {
  bool feasible;
  int64_t score;
  std::vector<int> codons;
};

// M3: greedy G1 (pure best codon per residue) or G2 (GC-steered tie-break).
GreedyResult greedy(const daa::dp::DpInstance& inst,
                    const daa::scores::ScoreTable& st,
                    const daa::automaton::AhoCorasick& ac,
                    const daa::motifs::MotifSet& ms,
                    bool gcSteered);

struct SampleResult {
  bool feasible;
  int64_t score;
  std::vector<int> codons;
  long long trials;
  long long wasted;
};

// M4: seeded sample-and-reject with a fixed trial budget.
SampleResult sampleAndReject(const daa::dp::DpInstance& inst,
                             const daa::scores::ScoreTable& st,
                             const daa::motifs::MotifSet& ms,
                             std::uint64_t seed,
                             long long budget,
                             double temperature);

}  // namespace baselines
}  // namespace daa
