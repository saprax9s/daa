cmake_minimum_required(VERSION 3.16)
project(codon_core VERSION 1.0.0 LANGUAGES CXX)

add_library(codon_core STATIC
  src/genetic_code.cpp
  src/prng.cpp
  src/scores.cpp
  src/motifs.cpp
  src/automaton.cpp
  src/dp.cpp
  src/baselines.cpp
  bench/bench.cpp
)

add_library(codon_core::codon_core ALIAS codon_core)

set_target_properties(codon_core PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)

target_include_directories(codon_core PUBLIC
  ${CMAKE_CURRENT_SOURCE_DIR}/src
  ${CMAKE_CURRENT_SOURCE_DIR}/bench
)

add_executable(codon_tests
  tests/selftest.cpp
)

target_link_libraries(codon_tests PRIVATE codon_core)

add_executable(codon_tool
  src/main.cpp
  tests/selftest.cpp
)

target_compile_definitions(codon_tool PRIVATE CODON_SELFTEST_NO_MAIN)

target_link_libraries(codon_tool PRIVATE codon_core)

set_target_properties(codon_tool PROPERTIES
  CXX_STANDARD 20
  CXX_STANDARD_REQUIRED ON
  CXX_EXTENSIONS OFF
)

enable_testing()

add_test(NAME selftest COMMAND codon_tests)

install(TARGETS codon_core codon_tests codon_tool
  DESTINATION lib
)
