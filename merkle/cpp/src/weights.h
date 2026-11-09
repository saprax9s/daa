#pragma once

#include <cstdint>
#include <vector>

namespace merkle {

// Weight model:
//   W_new(a) = (3 * touches(a) * 1000000 + 7 * W_old(a)) / 10   (integer division)
// updated once per epoch from touches only.

uint64_t updateWeight(uint64_t oldWeight, uint64_t touches);

// Hot set = top K scores, ties broken by lowest account index.
// Returns account indices (ascending) of the selected hot accounts.
std::vector<uint64_t> computeHotSet(const std::vector<uint64_t>& scores, uint64_t K);

}  // namespace merkle
