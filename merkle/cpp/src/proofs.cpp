#include "proofs.h"

#include <cstdio>
#include <map>
#include <string>

namespace merkle {

namespace {

Digest combine(const Digest& node, const ProofStep& s) {
    if (s.bit == 0) {
        // node is the left child -> parent = H(node, sibling)
        return hashInternal(node, s.siblingHash);
    } else {
        // node is the right child -> parent = H(sibling, node)
        return hashInternal(s.siblingHash, node);
    }
}

}  // namespace

MembershipProof generateMembership(const MerkleTree& t, int leafIdx) {
    MembershipProof p;
    p.leafIdx = leafIdx;
    p.leafHash = t.hashes[leafIdx];
    p.depth = t.depths[leafIdx];
    int cur = leafIdx;
    while (cur != t.root && cur >= 0) {
        int pId = t.nodes[cur].parent;
        if (pId < 0) break;
        bool isLeft = (t.nodes[pId].left == cur);
        int sibId = isLeft ? t.nodes[pId].right : t.nodes[pId].left;
        ProofStep s;
        s.bit = isLeft ? 0 : 1;
        s.siblingHash = t.hashes[sibId];
        p.steps.push_back(s);
        cur = pId;
    }
    return p;
}

bool verifyMembership(const MerkleTree& t, const MembershipProof& p) {
    Digest node = p.leafHash;
    for (const auto& s : p.steps) {
        node = combine(node, s);
    }
    return node == t.rootHash;
}

NonInclusionProof generateNonInclusion(const MerkleTree& t, int leafIdx) {
    // Adjacent-pair gap proof: full path from the lower adjacent leaf to the root,
    // plus the adjacency claim (upper leaf is lower+1 in sorted order).
    NonInclusionProof p;
    p.lowerIdx = leafIdx;
    p.upperIdx = leafIdx + 1;
    p.leafHash = t.hashes[leafIdx];
    p.depth = t.depths[leafIdx];
    int cur = leafIdx;
    while (cur != t.root && cur >= 0) {
        int pId = t.nodes[cur].parent;
        if (pId < 0) break;
        bool isLeft = (t.nodes[pId].left == cur);
        int sibId = isLeft ? t.nodes[pId].right : t.nodes[pId].left;
        ProofStep s;
        s.bit = isLeft ? 0 : 1;
        s.siblingHash = t.hashes[sibId];
        p.steps.push_back(s);
        cur = pId;
    }
    return p;
}

bool verifyNonInclusion(const MerkleTree& t, const NonInclusionProof& p) {
    if (p.upperIdx != p.lowerIdx + 1) return false;  // adjacency check
    Digest node = p.leafHash;
    for (const auto& s : p.steps) {
        node = combine(node, s);
    }
    return node == t.rootHash;
}

Multiproof generateMultiproof(const MerkleTree& t, const std::vector<int>& leafIds) {
    Multiproof mp;
    mp.rootHash = t.rootHash;
    std::map<std::string, uint64_t> idByHash;
    for (int leafId : leafIds) {
        MembershipProof p = generateMembership(t, leafId);
        std::vector<std::pair<uint64_t, uint8_t>> path;
        Digest node = p.leafHash;
        for (const auto& s : p.steps) {
            std::string key;
            for (int b = 0; b < 32; ++b) {
                char hex[3];
                std::snprintf(hex, 3, "%02x", s.siblingHash[b]);
                key += hex;
            }
            auto it = idByHash.find(key);
            uint64_t id;
            if (it == idByHash.end()) {
                id = (uint64_t)mp.sharedPool.size();
                idByHash[key] = id;
                mp.sharedPool.push_back(s.siblingHash);
            } else {
                id = it->second;
            }
            path.push_back({id, s.bit});
            node = combine(node, s);
        }
        (void)node;
        mp.leafPaths.push_back(path);
    }
    mp.uniqueHashCount = (size_t)mp.sharedPool.size();
    return mp;
}

bool verifyMultiproof(const Multiproof& mp, const MerkleTree& t, const std::vector<int>& leafIds) {
    if (mp.leafPaths.size() != leafIds.size()) return false;
    for (size_t i = 0; i < leafIds.size(); ++i) {
        Digest node = t.hashes[leafIds[i]];
        for (const auto& pr : mp.leafPaths[i]) {
            if (pr.first >= mp.sharedPool.size()) return false;
            node = combine(node, ProofStep{pr.second, mp.sharedPool[pr.first]});
        }
        if (node != t.rootHash) return false;
    }
    return true;
}

ProofSizeInfo proofSizeInfo(uint64_t hashes, uint64_t directionBits, uint64_t leafDataBytes) {
    ProofSizeInfo info;
    info.hashes = hashes;
    info.hashBytes = hashes * 32;
    info.bitBytes = (directionBits + 7) / 8;
    info.leafDataBytes = leafDataBytes;
    info.totalBytes = info.hashBytes + info.bitBytes + info.leafDataBytes;
    return info;
}

}  // namespace merkle
