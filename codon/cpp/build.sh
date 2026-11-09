#!/bin/bash
# Stage 1 build script for codon_core (Linux / macOS).
# If CMake is available, use the standard CMake flow; otherwise fall back to direct g++.

if command -v cmake >/dev/null 2>&1; then
  echo "Using CMake..."
  cmake -S . -B build
  cmake --build build --config Release
else
  echo "CMake not found, using direct g++..."
  mkdir -p build
  echo "Building codon_selftest..."
  g++ -std=c++20 -O2 -Wall -Wextra -Isrc -Itests -Ibench src/genetic_code.cpp src/prng.cpp src/scores.cpp src/motifs.cpp src/automaton.cpp src/dp.cpp src/baselines.cpp bench/bench.cpp tests/selftest.cpp -o build/codon_selftest
  echo "Building codon_tool..."
  g++ -std=c++20 -O2 -Wall -Wextra -Isrc -Itests -Ibench -DCODON_SELFTEST_NO_MAIN src/main.cpp src/genetic_code.cpp src/prng.cpp src/scores.cpp src/motifs.cpp src/automaton.cpp src/dp.cpp src/baselines.cpp bench/bench.cpp tests/selftest.cpp -o build/codon_tool
fi

echo "Build finished."
