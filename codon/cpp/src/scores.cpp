#include "scores.h"
#include "prng.h"

namespace daa {
namespace scores {

ScoreTable makeScoreTable(std::uint64_t seed) {
  ScoreTable st;
  st.codonScore.assign(64, 0);
  st.pairScore.assign(64, std::vector<int64_t>(64, 0));
  prng::SplitMix64 rng(seed);
  for (int i = 0; i < 64; ++i) {
    const std::uint64_t h = rng.nextUInt64();
    st.codonScore[i] = static_cast<int64_t>(static_cast<int64_t>(h % 2001) - 1000);
  }
  for (int i = 0; i < 64; ++i) {
    for (int j = 0; j < 64; ++j) {
      const std::uint64_t h = rng.nextUInt64();
      st.pairScore[i][j] = static_cast<int64_t>(static_cast<int64_t>(h % 1001) - 500);
    }
  }
  return st;
}

int64_t sequenceScore(const std::vector<int>& codons, const ScoreTable& st) {
  if (codons.empty()) return 0;
  int64_t total = 0;
  for (int c : codons) total += st.codonScore[c];
  for (std::size_t i = 1; i < codons.size(); ++i) {
    total += st.pairScore[codons[i - 1]][codons[i]];
  }
  return total;
}

}  // namespace scores
}  // namespace daa
