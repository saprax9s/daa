#include "sha256.h"

namespace merkle {

namespace {
inline uint32_t rotr32(uint32_t x, int n) { return (x >> n) | (x << (32 - n)); }
}  // namespace

Sha256::Sha256() { reset(); }

void Sha256::reset() {
    state_[0] = 0x6a09e667u; state_[1] = 0xbb67ae85u; state_[2] = 0x3c6ef372u;
    state_[3] = 0xa54ff53au; state_[4] = 0x510e527fu; state_[5] = 0x9b05688cu;
    state_[6] = 0x1f83d9abu; state_[7] = 0x5be0cd19u;
    bitLen_ = 0;
    bufferLen_ = 0;
}

void Sha256::update(const uint8_t* data, size_t len) {
    bitLen_ += static_cast<uint64_t>(len) * 8u;
    while (len > 0) {
        size_t take = 64 - bufferLen_;
        if (take > len) take = len;
        for (size_t i = 0; i < take; ++i) buffer_[bufferLen_ + i] = data[i];
        bufferLen_ += take;
        data += take;
        len -= take;
        if (bufferLen_ == 64) {
            transform(buffer_);
            bufferLen_ = 0;
        }
    }
}

void Sha256::update(const std::string& s) { update(reinterpret_cast<const uint8_t*>(s.data()), s.size()); }

void Sha256::final(uint8_t out[32]) {
    uint64_t bits = bitLen_;
    size_t len = bufferLen_;
    buffer_[len++] = 0x80;
    if (len > 56) {
        while (len < 64) buffer_[len++] = 0x00;
        transform(buffer_);
        len = 0;
    }
    while (len < 56) buffer_[len++] = 0x00;
    for (int i = 0; i < 8; ++i) {
        buffer_[56 + i] = static_cast<uint8_t>((bits >> (56 - 8 * i)) & 0xFFu);
    }
    transform(buffer_);
    for (int i = 0; i < 8; ++i) {
        out[4 * i] = static_cast<uint8_t>(state_[i] >> 24);
        out[4 * i + 1] = static_cast<uint8_t>(state_[i] >> 16);
        out[4 * i + 2] = static_cast<uint8_t>(state_[i] >> 8);
        out[4 * i + 3] = static_cast<uint8_t>(state_[i]);
    }
}

void Sha256::transform(const uint8_t* p) {
    static const uint32_t K[64] = {
        0x428a2f98u, 0x71374491u, 0xb5c0fbcfu, 0xe9b5dba5u, 0x3956c25bu, 0x59f111f1u, 0x923f82a4u, 0xab1c5ed5u,
        0xd807aa98u, 0x12835b01u, 0x243185beu, 0x550c7dc3u, 0x72be5d74u, 0x80deb1feu, 0x9bdc06a7u, 0xc19bf174u,
        0xe49b69c1u, 0xefbe4786u, 0x0fc19dc6u, 0x240ca1ccu, 0x2de92c6fu, 0x4a7484aau, 0x5cb0a9dcu, 0x76f988dau,
        0x983e5152u, 0xa831c66du, 0xb00327c8u, 0xbf597fc7u, 0xc6e00bf3u, 0xd5a79147u, 0x06ca6351u, 0x14292967u,
        0x27b70a85u, 0x2e1b2138u, 0x4d2c6dfcu, 0x53380d13u, 0x650a7354u, 0x766a0abbu, 0x81c2c92eu, 0x92722c85u,
        0xa2bfe8a1u, 0xa81a664bu, 0xc24b8b70u, 0xc76c51a3u, 0xd192e819u, 0xd6990624u, 0xf40e3585u, 0x106aa070u,
        0x19a4c116u, 0x1e376c08u, 0x2748774cu, 0x34b0bcb5u, 0x391c0cb3u, 0x4ed8aa4au, 0x5b9cca4fu, 0x682e6ff3u,
        0x748f82eeu, 0x78a5636fu, 0x84c87814u, 0x8cc70208u, 0x90befffau, 0xa4506cebu, 0xbef9a3f7u, 0xc67178f2u};
    uint32_t w[64];
    for (int i = 0; i < 16; ++i) {
        w[i] = (static_cast<uint32_t>(p[4 * i]) << 24) | (static_cast<uint32_t>(p[4 * i + 1]) << 16) |
               (static_cast<uint32_t>(p[4 * i + 2]) << 8) | static_cast<uint32_t>(p[4 * i + 3]);
    }
    for (int i = 16; i < 64; ++i) {
        uint32_t s0 = rotr32(w[i - 15], 7) ^ rotr32(w[i - 15], 18) ^ (w[i - 15] >> 3);
        uint32_t s1 = rotr32(w[i - 2], 17) ^ rotr32(w[i - 2], 19) ^ (w[i - 2] >> 10);
        w[i] = w[i - 16] + s0 + w[i - 7] + s1;
    }
    uint32_t a = state_[0], b = state_[1], c = state_[2], d = state_[3];
    uint32_t e = state_[4], f = state_[5], g = state_[6], h = state_[7];
    for (int i = 0; i < 64; ++i) {
        uint32_t S1 = rotr32(e, 6) ^ rotr32(e, 11) ^ rotr32(e, 25);
        uint32_t ch = (e & f) ^ ((~e) & g);
        uint32_t t1 = h + S1 + ch + K[i] + w[i];
        uint32_t S0 = rotr32(a, 2) ^ rotr32(a, 13) ^ rotr32(a, 22);
        uint32_t maj = (a & b) ^ (a & c) ^ (b & c);
        uint32_t t2 = S0 + maj;
        h = g; g = f; f = e; e = d + t1; d = c; c = b; b = a; a = t1 + t2;
    }
    state_[0] += a; state_[1] += b; state_[2] += c; state_[3] += d;
    state_[4] += e; state_[5] += f; state_[6] += g; state_[7] += h;
}

std::array<uint8_t, 32> sha256Bytes(const uint8_t* data, size_t len) {
    Sha256 ctx;
    ctx.update(data, len);
    std::array<uint8_t, 32> out{};
    ctx.final(out.data());
    return out;
}

std::array<uint8_t, 32> sha256Bytes(const std::string& s) { return sha256Bytes(reinterpret_cast<const uint8_t*>(s.data()), s.size()); }

std::string sha256Hex(const std::string& s) {
    auto d = sha256Bytes(s);
    static const char* hexd = "0123456789abcdef";
    std::string out;
    out.reserve(64);
    for (int i = 0; i < 32; ++i) {
        out.push_back(hexd[d[i] >> 4]);
        out.push_back(hexd[d[i] & 0xf]);
    }
    return out;
}

}  // namespace merkle
