cmake_minimum_required(VERSION 3.16)
project(merkle VERSION 1.0.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

set(MERKLE_CORE_SOURCES
    src/sha256.cpp
    src/hashing.cpp
    src/weights.cpp
    src/trees.cpp
    src/proofs.cpp
    src/workloads.cpp
    src/pipeline.cpp
    bench/bench.cpp
    tests/selftest.cpp
)

add_library(merkle_core STATIC ${MERKLE_CORE_SOURCES})
target_include_directories(merkle_core PUBLIC "${CMAKE_SOURCE_DIR}/src" "${CMAKE_SOURCE_DIR}/tests" "${CMAKE_SOURCE_DIR}/bench")

add_executable(merkle_tool src/main.cpp)
target_link_libraries(merkle_tool PRIVATE merkle_core)

add_executable(selftest tests/selftest.cpp)
target_compile_definitions(selftest PRIVATE SELFTEST_STANDALONE_MAIN)
target_link_libraries(selftest PRIVATE merkle_core)

enable_testing()
add_test(NAME selftest COMMAND selftest)
