#pragma once

#include <cstdint>
#include <vector>

#include "automaton.h"
#include "genetic_code.h"
#include "motifs.h"
#include "scores.h"

namespace daa {
namespace dp {

struct DpInstance {
  int n = 0;  // protein length
  std::vector<int> protein;                // amino acid indices
  std::vector<std::vector<int>> allowed;  // per-residue allowed codon indices
  int64_t lo = 0;  // global GC lower bound
  int64_t hi = 0;  // global GC upper bound
  int windowK = 0;  // window GC enabled with window size k (0 = off)
  int windowLo = 0;
  int windowHi = 0;
};

struct DpFlags {
  bool gcPruning = false;  // D1: GC-bound pruning
  bool dominance = false;  // D2: D1 + dominance pruning (only valid when lo <= 0)
};

struct DpResult {
  bool feasible = false;
  int64_t score = 0;
  std::vector<int> codons;  // chosen codon indices (size n); empty if infeasible
  uint64_t totalStates = 0;
  uint64_t peakStates = 0;
};

// Exact DP over layers: state = (last codon, automaton state, window GC vector, GC count).
DpResult solve(const daa::automaton::AhoCorasick& ac,
              const daa::scores::ScoreTable& st,
              const daa::motifs::MotifSet& ms,
              const daa::dp::DpInstance& inst,
              const daa::dp::DpFlags& flags);

}  // namespace dp
}  // namespace daa
