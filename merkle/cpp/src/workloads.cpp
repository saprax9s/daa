#include "workloads.h"

#include <algorithm>
#include <cmath>

namespace merkle {

uint64_t splitmix64(uint64_t& state) {
    state += 0x9E3779B97F4A7C15ULL;
    uint64_t z = state;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    return z ^ (z >> 31);
}

void splitmix64Seed(uint64_t& state, uint64_t seed) { state = seed; }

std::vector<std::vector<uint64_t>> generateEpochTouches(uint64_t seed, uint64_t N, double s, uint64_t blocksPerEpoch, uint64_t txPerBlock) {
    // Zipf access model: P(rank) = (rank+1)^(-s) / H_N(s), rank = 0-based account index.
    std::vector<double> cum(N, 0.0);
    double total = 0.0;
    std::vector<double> p(N, 0.0);
    for (uint64_t i = 0; i < N; ++i) {
        p[i] = std::pow((double)(i + 1), -s);
        total += p[i];
    }
    double acc = 0.0;
    for (uint64_t i = 0; i < N; ++i) {
        acc += p[i] / total;
        cum[i] = acc;
    }
    std::vector<std::vector<uint64_t>> epochs;
    epochs.reserve(blocksPerEpoch);
    uint64_t state = seed;
    for (uint64_t b = 0; b < blocksPerEpoch; ++b) {
        std::vector<uint64_t> touches(N, 0);
        for (uint64_t t = 0; t < txPerBlock; ++t) {
            double u1 = (double)(splitmix64(state) % 1000000) / 1000000.0;
            double u2 = (double)(splitmix64(state) % 1000000) / 1000000.0;
            uint64_t a = (uint64_t)(std::lower_bound(cum.begin(), cum.end(), u1) - cum.begin());
            uint64_t c = (uint64_t)(std::lower_bound(cum.begin(), cum.end(), u2) - cum.begin());
            touches[a] += 1;
            touches[c] += 1;
        }
        epochs.push_back(touches);
    }
    return epochs;
}

std::vector<uint64_t> applyCascade(const std::vector<uint64_t>& weights, uint64_t seed, int mode) {
    // Attacker controls the first attackerCount leaves (already sorted by key).
    // mode 0: geometric growth of touches; mode 1: Fibonacci additions.
    std::vector<uint64_t> w = weights;
    int cnt = (int)std::min((size_t)64, (size_t)w.size());
    uint64_t state = seed;
    if (mode == 0) {
        // geometric: touches_k = base_k * (1 + growth^k)
        double growth = 1.0 + (double)(splitmix64(state) % 10) / 10.0;
        for (int i = 0; i < cnt; ++i) {
            uint64_t base = (i < (int)w.size()) ? w[i] : 1;
            uint64_t add = (uint64_t)((double)base * growth * (double)(splitmix64(state) % 4 + 1) + 1.0);
            w[i] = base + add;
        }
    } else {
        // Fibonacci cascade: touches_i += fib(i)
        uint64_t a = 1, b = 1;
        for (int i = 0; i < cnt; ++i) {
            uint64_t fib = (i < 2) ? 1 : a + b;
            a = b;
            b = fib;
            if (i < (int)w.size()) w[i] = w[i] + fib;
        }
    }
    return w;
}

std::vector<uint64_t> attackerBudgetTouches(const std::vector<uint64_t>& hotWeights, uint64_t seed, uint64_t perBlockCap, uint64_t perEpochTotal) {
    // Budget-limited attacker: touches are spread over the hot set under a
    // per-block cap and per-epoch total budget (gas-derived).
    std::vector<uint64_t> touches(hotWeights.size(), 0);
    uint64_t state = seed;
    uint64_t spent = 0;
    for (uint64_t i = 0; i < touches.size() && spent < perEpochTotal; ++i) {
        uint64_t u = splitmix64(state) % 100;
        uint64_t t = (u < 60) ? 0 : (1 + (splitmix64(state) % (perBlockCap + 1)));
        if (spent + t > perEpochTotal) t = perEpochTotal - spent;
        touches[i] = t;
        spent += t;
    }
    return touches;
}

std::vector<uint64_t> applyGriefing(const std::vector<uint64_t>& weights, uint64_t seed, double intensity) {
    std::vector<uint64_t> w = weights;
    int n = (int)w.size();
    uint64_t state = seed;
    for (int i = 0; i < n; ++i) {
        if (i == 0) continue;  // skip the largest
        uint64_t add = (uint64_t)(1 + (double)(splitmix64(state) % 100) / 100.0 * intensity * 10.0);
        w[i] = w[i] + add;
    }
    return w;
}

}  // namespace merkle
