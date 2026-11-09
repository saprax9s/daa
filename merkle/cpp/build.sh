#!/usr/bin/env bash
set -e
mkdir -p build
if command -v cmake >/dev/null 2>&1; then
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --config Release
else
    g++ -std=c++20 -O2 -Wall -Wextra -Isrc -Itests -Ibench src/sha256.cpp src/hashing.cpp src/weights.cpp src/trees.cpp src/proofs.cpp src/workloads.cpp src/pipeline.cpp bench/bench.cpp tests/selftest.cpp src/main.cpp -o build/merkle_tool
fi
echo "Build finished."
