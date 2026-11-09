#include "hashing.h"

#include <cstdio>
#include "sha256.h"

namespace merkle {

namespace {
std::string buildMessage(const uint8_t tag, const std::string& body) {
    std::string msg;
    msg.push_back(static_cast<char>(tag));
    msg += body;
    return msg;
}
}  // namespace

Digest hashAccount(uint64_t key, uint64_t weight) {
    std::string body;
    for (int i = 0; i < 8; ++i) body.push_back(static_cast<char>((key >> (8 * i)) & 0xFFu));
    for (int i = 0; i < 8; ++i) body.push_back(static_cast<char>((weight >> (8 * i)) & 0xFFu));
    auto d = sha256Bytes(buildMessage(static_cast<uint8_t>(0x00), body));
    return d;
}

Digest hashInternal(const Digest& left, const Digest& right) {
    std::string msg;
    msg.push_back(static_cast<char>(0x01));
    msg.append(reinterpret_cast<const char*>(left.data()), 32);
    msg.append(reinterpret_cast<const char*>(right.data()), 32);
    auto d = sha256Bytes(msg);
    return d;
}

Digest hashCold(uint64_t bucketWeight, uint64_t bucketLeafCount) {
    std::string body;
    for (int i = 0; i < 8; ++i) body.push_back(static_cast<char>((bucketWeight >> (8 * i)) & 0xFFu));
    for (int i = 0; i < 8; ++i) body.push_back(static_cast<char>((bucketLeafCount >> (8 * i)) & 0xFFu));
    auto d = sha256Bytes(buildMessage(static_cast<uint8_t>(0x02), body));
    return d;
}

Digest hashTop(const Digest& rootDigest, uint64_t totalWeight, uint64_t hotCount) {
    std::string body;
    body.append(reinterpret_cast<const char*>(rootDigest.data()), 32);
    for (int i = 0; i < 8; ++i) body.push_back(static_cast<char>((totalWeight >> (8 * i)) & 0xFFu));
    for (int i = 0; i < 8; ++i) body.push_back(static_cast<char>((hotCount >> (8 * i)) & 0xFFu));
    auto d = sha256Bytes(buildMessage(static_cast<uint8_t>(0x03), body));
    return d;
}

std::string digestHex(const Digest& d) {
    char buf[65];
    for (int i = 0; i < 32; ++i) {
        snprintf(buf + 2 * i, 3, "%02x", d[i]);
    }
    return std::string(buf, 64);
}

std::string digestToHex(const uint8_t* d) {
    Digest dg{};
    for (int i = 0; i < 32; ++i) dg[i] = d[i];
    return digestHex(dg);
}

}  // namespace merkle
