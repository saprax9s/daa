#pragma once

#include <array>
#include <cstdint>
#include <string>

namespace merkle {

class Sha256 {
public:
    Sha256();
    void reset();
    void update(const uint8_t* data, size_t len);
    void update(const std::string& s);
    void final(uint8_t out[32]);

private:
    void transform(const uint8_t* p);

    uint32_t state_[8];
    uint64_t bitLen_;
    uint8_t buffer_[64];
    size_t bufferLen_;
};

std::array<uint8_t, 32> sha256Bytes(const uint8_t* data, size_t len);
std::array<uint8_t, 32> sha256Bytes(const std::string& s);
std::string sha256Hex(const std::string& s);

}  // namespace merkle
