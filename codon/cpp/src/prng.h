#pragma once

#include <cstdint>

namespace daa {
namespace prng {

// Deterministic splitmix64 pseudo-random generator.
class SplitMix64 {
 public:
  explicit SplitMix64(std::uint64_t seed) : state_(seed) {}

  std::uint64_t nextUInt64() {
    state_ += 0x9E3779B97F4A7C15ULL;
    std::uint64_t z = state_;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
  }

  // Uniform integer in [0, n).
  int below(std::uint64_t n) { return static_cast<int>(nextUInt64() % n); }

  // Uniform double in [0, 1).
  double nextDouble() {
    return (nextUInt64() >> 11) * (1.0 / 9007199254740992.0);
  }

 private:
  std::uint64_t state_;
};

}  // namespace prng
}  // namespace daa
