#include <cstdio>
#include <string>
#include <vector>

#include "bench.h"
#include "selftest.h"
#include "trees.h"

int main(int argc, char** argv) {
    for (int i = 1; i < argc; ++i) {
        std::string arg = argv[i];
        if (arg == "--selftest") {
            merkle::SelfTestResult r = merkle::runSelfTest();
            std::printf("selftest: %s (%d passed, %d failed)\n", r.details.c_str(), r.checksPassed, r.checksFailed);
            return (r.checksFailed == 0) ? 0 : 1;
        }
        if (arg == "--bench") {
            std::string outFile = "../results/merkle_results.csv";
            if (i + 1 < argc && argv[i + 1][0] != '-') {
                outFile = argv[++i];
            }
            return merkle::runExperiments(outFile);
        }
        if (arg == "--cost") {
            // usage: merkle_tool --cost <leaves> <seed>
            if (i + 2 >= argc) {
                std::printf("usage: merkle_tool --selftest | --bench [csv_path] | --cost <leaves> <seed>\n");
                return 2;
            }
            int n = std::atoi(argv[++i]);
            long seed = std::atol(argv[++i]);
            std::vector<uint64_t> w((size_t)n, 0);
            uint64_t state = (uint64_t)seed;
            for (int j = 0; j < n; ++j) {
                state += 0x9E3779B97F4A7C15ULL;
                uint64_t z = state;
                z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
                z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
                z = z ^ (z >> 31);
                w[(size_t)j] = z % 100ULL + 1ULL;
            }
            merkle::TreeResult r = merkle::solveOptimal(w);
            std::printf("optimal: leaves=%d cost=%llu maxDepth=%llu\n", n, (unsigned long long)r.cost, (unsigned long long)r.maxDepth);
            return 0;
        }
    }
    std::printf("usage: merkle_tool --selftest | --bench [csv_path] | --cost <leaves> <seed>\n");
    return 2;
}
