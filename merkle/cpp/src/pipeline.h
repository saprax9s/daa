#pragma once

#include <cstdint>
#include <vector>

namespace merkle {

// Tarjan SCC for wash-trade detection
struct DirectedGraph {
    int numVertices = 0;
    std::vector<std::vector<int>> adj;
};

struct SccResult {
    std::vector<std::vector<int>> components;
    int cycleCount = 0; // number of components with size > 1 or self-loops
    bool hasWashTradeCycle = false;
};

SccResult runTarjanSCC(const DirectedGraph& g);
bool bruteForceHasCycle(const DirectedGraph& g);

// Pre-block Knapsack DP vs Greedy fee-per-gas
struct Tx {
    int id;
    int account;
    int nonce;
    uint64_t gas;
    uint64_t fee;
};

struct BlockSelectionResult {
    uint64_t totalFee = 0;
    uint64_t totalGas = 0;
    std::vector<int> selectedTxIds;
};

// Exact DP knapsack respecting nonce-prefix constraints per account
BlockSelectionResult solveBlockKnapsackDP(const std::vector<Tx>& txs, uint64_t gasLimit);

// Greedy heuristic: fee-per-gas picking available nonces
BlockSelectionResult solveBlockGreedy(const std::vector<Tx>& txs, uint64_t gasLimit);

}  // namespace merkle

