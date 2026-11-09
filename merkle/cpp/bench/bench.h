#pragma once

#include <string>

namespace merkle {

// Runs all experiments E1 to E9, outputs results to results/merkle_results.csv
int runExperiments(const std::string& outputFile);

}  // namespace merkle

