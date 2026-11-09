@echo off
setlocal
if not exist build mkdir build

where cmake >nul 2>nul
if %errorlevel%==0 (
    echo Using CMake...
    cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
    cmake --build build --config Release
) else (
    echo CMake not found, using direct g++...
    g++ -std=c++20 -O2 -Wall -Wextra -Isrc -Itests -Ibench src/sha256.cpp src/hashing.cpp src/weights.cpp src/trees.cpp src/proofs.cpp src/workloads.cpp src/pipeline.cpp bench/bench.cpp tests/selftest.cpp src/main.cpp -o build\merkle_tool.exe
)
echo Build finished.
