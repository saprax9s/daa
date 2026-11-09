#pragma once

#include <cstdint>
#include <vector>

#include "hashing.h"

namespace merkle {

struct TreeResult {
    uint64_t cost = 0;      // weighted external path length: sum_i w_i * depth_i
    uint64_t maxDepth = 0;  // maximum leaf depth (root = depth 0)
    uint64_t leafCount = 0;
};

// B1/B2: balanced binary tree over sorted leaves.
TreeResult solveBalanced(const std::vector<uint64_t>& w);

// B4: Garsia-Wachs optimal alphabetic tree, simple O(K^2) variant
//     (repeatedly merge the adjacent pair with the minimal weight sum).
TreeResult solveOptimal(const std::vector<uint64_t>& w);

// B5: Knuth interval DP with monotone split-point optimization, O(K^2).
TreeResult solveKnuth(const std::vector<uint64_t>& w);

// B6: depth-capped alphabetic DP DP[h][i][j], O(K^2 * Dmax).
TreeResult solveDepthCapped(const std::vector<uint64_t>& w, uint64_t Dmax);

// B3: Mehlhorn simple version (greedy adjacent-pair merge), O(K^2).
TreeResult solveMehlhorn(const std::vector<uint64_t>& w);

// B7: weight-clamped Garsia-Wachs: raise every weight to at least W_total / 2^c.
TreeResult solveClamped(const std::vector<uint64_t>& w, int c);

// B8: Huffman-layout hot tier + balanced cold tier (order NOT preserved).
TreeResult solveHuffman(const std::vector<uint64_t>& hotWeights, const std::vector<uint64_t>& coldBucketWeights);

// B9: static top-L caching, analytic: cost = sum_i w_i * max(0, depth_i - L).
uint64_t cachedCost(const std::vector<uint64_t>& w, uint64_t L, const std::vector<uint64_t>& depths);

// Proof-supporting tree build (optimal greedy construction) over sorted weights.
struct TreeNode {
    bool isLeaf;
    int leafIdx;
    int left;
    int right;
    int parent;
};

struct MerkleTree {
    std::vector<TreeNode> nodes;
    int root = -1;
    std::vector<uint64_t> depths;  // per leaf
    std::vector<Digest> hashes;    // per node
    Digest rootHash{};
};

MerkleTree buildOptimalTree(const std::vector<uint64_t>& w);

}  // namespace merkle
