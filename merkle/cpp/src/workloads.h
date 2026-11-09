#pragma once

#include <cstdint>
#include <vector>

namespace merkle {

// Deterministic splitmix64 PRNG.
uint64_t splitmix64(uint64_t& state);

void splitmix64Seed(uint64_t& state, uint64_t seed);

// Zipf account access generator: exponent s, blocks of T transactions each
// touching 2 accounts; epochs of E blocks. Returns per-epoch touch vectors
// (index = account id, value = touches in that epoch).
std::vector<std::vector<uint64_t>> generateEpochTouches(uint64_t seed, uint64_t N, double s, uint64_t blocksPerEpoch, uint64_t txPerBlock);

struct AdversaryConfig {
    int mode = 0;  // 0 = geometric/Fibonacci weight cascade, 1 = budget-limited, 2 = griefing
    uint64_t attackerCount = 0;
};

// (a) geometric / Fibonacci weight cascade on a set of attacker-controlled leaves.
std::vector<uint64_t> applyCascade(const std::vector<uint64_t>& weights, uint64_t seed, int mode);

// (b) budget-limited attacker: per-block touch cap + per-epoch total budget.
std::vector<uint64_t> attackerBudgetTouches(const std::vector<uint64_t>& hotWeights, uint64_t seed, uint64_t perBlockCap, uint64_t perEpochTotal);

// (c) griefing: push honest hot accounts deeper by adding touches to them.
std::vector<uint64_t> applyGriefing(const std::vector<uint64_t>& weights, uint64_t seed, double intensity);

}  // namespace merkle
