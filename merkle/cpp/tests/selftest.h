#pragma once

#include <string>

namespace merkle {

struct SelfTestResult {
    int checksPassed = 0;
    int checksFailed = 0;
    std::string details;
};

SelfTestResult runSelfTest();

}  // namespace merkle
