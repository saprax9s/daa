#include "weights.h"

#include <algorithm>

namespace merkle {

uint64_t updateWeight(uint64_t oldWeight, uint64_t touches) {
    return (3ULL * touches * 1000000ULL + 7ULL * oldWeight) / 10ULL;
}

std::vector<uint64_t> computeHotSet(const std::vector<uint64_t>& scores, uint64_t K) {
    // entries: (score desc, index asc)
    std::vector<std::pair<uint64_t, uint64_t>> entries;
    entries.reserve(scores.size());
    for (uint64_t i = 0; i < scores.size(); ++i) {
        entries.push_back({scores[i], i});
    }
    std::sort(entries.begin(), entries.end(), [](const std::pair<uint64_t, uint64_t>& a, const std::pair<uint64_t, uint64_t>& b) {
        if (a.first != b.first) return a.first > b.first;  // score descending
        return a.second < b.second;                        // ties -> lowest index first
    });
    std::vector<uint64_t> hot;
    if (K > entries.size()) K = entries.size();
    for (uint64_t i = 0; i < K; ++i) hot.push_back(entries[i].second);
    std::sort(hot.begin(), hot.end());
    return hot;
}

}  // namespace merkle
