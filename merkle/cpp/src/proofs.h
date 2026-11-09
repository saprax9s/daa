#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "hashing.h"
#include "trees.h"

namespace merkle {

struct ProofStep {
    uint8_t bit = 0;  // 0 = current node is the left child, 1 = current node is the right child
    Digest siblingHash{};
};

struct MembershipProof {
    int leafIdx = -1;
    Digest leafHash{};
    uint64_t depth = 0;
    std::vector<ProofStep> steps;  // ordered from leaf (deepest) to root
};

MembershipProof generateMembership(const MerkleTree& t, int leafIdx);

bool verifyMembership(const MerkleTree& t, const MembershipProof& p);

struct NonInclusionProof {
    int lowerIdx = -1;  // adjacent lower leaf
    int upperIdx = -1;  // adjacent upper leaf, must be lowerIdx + 1
    Digest leafHash{};
    uint64_t depth = 0;
    std::vector<ProofStep> steps;  // full path from lower leaf up to the root
};

NonInclusionProof generateNonInclusion(const MerkleTree& t, int leafIdx);

bool verifyNonInclusion(const MerkleTree& t, const NonInclusionProof& p);

struct Multiproof {
    Digest rootHash{};
    std::vector<std::vector<std::pair<uint64_t, uint8_t>>> leafPaths;  // per leaf: (sharedHashId, bit) in order
    std::vector<Digest> sharedPool;
    size_t uniqueHashCount = 0;
};

Multiproof generateMultiproof(const MerkleTree& t, const std::vector<int>& leafIds);

bool verifyMultiproof(const Multiproof& mp, const MerkleTree& t, const std::vector<int>& leafIds);

struct ProofSizeInfo {
    uint64_t hashes = 0;
    uint64_t hashBytes = 0;
    uint64_t bitBytes = 0;
    uint64_t leafDataBytes = 0;
    uint64_t totalBytes = 0;
};

ProofSizeInfo proofSizeInfo(uint64_t hashes, uint64_t directionBits, uint64_t leafDataBytes);

}  // namespace merkle
