#pragma once

#include <array>
#include <cstdint>
#include <string>

namespace merkle {

using Digest = std::array<uint8_t, 32>;

// Domain-separated hashing with a 1-byte tag (FIPS-style domain separation).
//
// Concrete scheme (documented, deterministic):
//   account leaf  (tag 0x00): SHA-256 over [0x00] ++ u64le(key) ++ u64le(weight)
//   internal node (tag 0x01): SHA-256 over [0x01] ++ leftDigest(32 raw bytes) ++ rightDigest(32 raw bytes)
//   cold bucket  (tag 0x02): SHA-256 over [0x02] ++ u64le(bucketWeight) ++ u64le(bucketLeafCount)
//   top commit   (tag 0x03): SHA-256 over [0x03] ++ rootDigest(32 raw bytes) ++ u64le(totalWeight) ++ u64le(hotCount)
//
// uint64 fields are encoded little-endian; digest bytes are used raw (big-endian byte order).

Digest hashAccount(uint64_t key, uint64_t weight);
Digest hashInternal(const Digest& left, const Digest& right);
Digest hashCold(uint64_t bucketWeight, uint64_t bucketLeafCount);
Digest hashTop(const Digest& rootDigest, uint64_t totalWeight, uint64_t hotCount);

std::string digestHex(const Digest& d);
std::string digestToHex(const uint8_t* d);

}  // namespace merkle
