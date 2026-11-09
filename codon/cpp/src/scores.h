#pragma once

#include <cstdint>
#include <vector>

namespace daa {
namespace scores {

struct ScoreTable {
  std::vector<int64_t> codonScore;               // size 64, per-codon additive score
  std::vector<std::vector<int64_t>> pairScore; // size 64 x 64, adjacent codon-pair scores
};

// Creates a seeded synthetic integer score table (deterministic).
ScoreTable makeScoreTable(std::uint64_t seed);

// Score of a codon sequence = sum of codon scores + sum of adjacent pair scores.
int64_t sequenceScore(const std::vector<int>& codons, const ScoreTable& st);

}  // namespace scores
}  // namespace daa
