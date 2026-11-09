#pragma once

#include <string>

namespace daa {

// Generates canonical test and benchmark instances to instances/*.json
int generateInstances(const std::string& instancesDir);

// Runs all experiments E1 to E8 and writes to results/codon_results.csv
int runCodonBenchmarks(const std::string& outputFile);

}  // namespace daa

