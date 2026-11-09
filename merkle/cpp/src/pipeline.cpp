#include "pipeline.h"

#include <algorithm>
#include <functional>
#include <stack>
#include <map>

namespace merkle {

SccResult runTarjanSCC(const DirectedGraph& g) {
    SccResult res;
    int n = g.numVertices;
    std::vector<int> dfn(n, -1);
    std::vector<int> low(n, -1);
    std::vector<bool> inStack(n, false);
    std::stack<int> st;
    int timer = 0;

    std::function<void(int)> dfs = [&](int u) {
        dfn[u] = low[u] = ++timer;
        st.push(u);
        inStack[u] = true;

        for (int v : g.adj[u]) {
            if (dfn[v] == -1) {
                dfs(v);
                low[u] = std::min(low[u], low[v]);
            } else if (inStack[v]) {
                low[u] = std::min(low[u], dfn[v]);
            }
        }

        if (low[u] == dfn[u]) {
            std::vector<int> comp;
            while (true) {
                int v = st.top();
                st.pop();
                inStack[v] = false;
                comp.push_back(v);
                if (u == v) break;
            }
            bool isCycle = false;
            if (comp.size() > 1) {
                isCycle = true;
            } else {
                // Check self loop
                for (int v : g.adj[u]) {
                    if (v == u) {
                        isCycle = true;
                        break;
                    }
                }
            }
            if (isCycle) {
                res.cycleCount++;
                res.hasWashTradeCycle = true;
            }
            res.components.push_back(comp);
        }
    };

    for (int i = 0; i < n; ++i) {
        if (dfn[i] == -1) dfs(i);
    }
    return res;
}

bool bruteForceHasCycle(const DirectedGraph& g) {
    int n = g.numVertices;
    std::vector<int> visited(n, 0); // 0 = unvisited, 1 = visiting, 2 = done

    std::function<bool(int)> dfs = [&](int u) -> bool {
        visited[u] = 1;
        for (int v : g.adj[u]) {
            if (visited[v] == 1) return true;
            if (visited[v] == 0 && dfs(v)) return true;
        }
        visited[u] = 2;
        return false;
    };

    for (int i = 0; i < n; ++i) {
        if (visited[i] == 0) {
            if (dfs(i)) return true;
        }
    }
    return false;
}

BlockSelectionResult solveBlockKnapsackDP(const std::vector<Tx>& txs, uint64_t gasLimit) {
    // Group transactions by account, sorted by nonce
    std::map<int, std::vector<Tx>> byAccount;
    for (const auto& tx : txs) {
        byAccount[tx.account].push_back(tx);
    }
    for (auto& pair : byAccount) {
        std::sort(pair.second.begin(), pair.second.end(), [](const Tx& a, const Tx& b) {
            return a.nonce < b.nonce;
        });
    }

    // For each account, create mutually exclusive prefixes (prefix length k: 1..len)
    // with cumulative gas and cumulative fee.
    struct PrefixOption {
        uint64_t gas;
        uint64_t fee;
        std::vector<int> txIds;
    };
    std::vector<std::vector<PrefixOption>> accountGroups;
    for (const auto& pair : byAccount) {
        std::vector<PrefixOption> options;
        uint64_t curGas = 0;
        uint64_t curFee = 0;
        std::vector<int> curIds;
        for (const auto& tx : pair.second) {
            curGas += tx.gas;
            curFee += tx.fee;
            curIds.push_back(tx.id);
            if (curGas <= gasLimit) {
                options.push_back(PrefixOption{curGas, curFee, curIds});
            }
        }
        if (!options.empty()) {
            accountGroups.push_back(options);
        }
    }

    // Multiple-choice knapsack DP: dp[g] = best fee using at most g gas
    std::vector<uint64_t> dp(gasLimit + 1, 0);
    // Track choices for reconstruction
    struct Choice {
        int groupIdx;
        int optIdx; // -1 if skipped
        uint64_t prevGas;
    };
    std::vector<std::vector<Choice>> parent(accountGroups.size(), std::vector<Choice>(gasLimit + 1, {-1, -1, 0}));

    for (size_t g = 0; g < accountGroups.size(); ++g) {
        std::vector<uint64_t> nextDp = dp;
        for (uint64_t w = 0; w <= gasLimit; ++w) {
            parent[g][w] = Choice{(int)g, -1, w};
        }
        for (size_t optIdx = 0; optIdx < accountGroups[g].size(); ++optIdx) {
            const auto& opt = accountGroups[g][optIdx];
            for (uint64_t w = opt.gas; w <= gasLimit; ++w) {
                uint64_t cand = dp[w - opt.gas] + opt.fee;
                if (cand > nextDp[w]) {
                    nextDp[w] = cand;
                    parent[g][w] = Choice{(int)g, (int)optIdx, w - opt.gas};
                }
            }
        }
        dp = std::move(nextDp);
    }

    uint64_t bestFee = 0;
    uint64_t bestGas = 0;
    for (uint64_t w = 0; w <= gasLimit; ++w) {
        if (dp[w] > bestFee) {
            bestFee = dp[w];
            bestGas = w;
        }
    }

    BlockSelectionResult res;
    res.totalFee = bestFee;
    res.totalGas = bestGas;
    if (accountGroups.empty()) return res;

    // Backtrack choices
    uint64_t curW = bestGas;
    for (int g = (int)accountGroups.size() - 1; g >= 0; --g) {
        Choice ch = parent[g][curW];
        if (ch.optIdx >= 0) {
            const auto& opt = accountGroups[g][ch.optIdx];
            res.selectedTxIds.insert(res.selectedTxIds.end(), opt.txIds.begin(), opt.txIds.end());
            curW = ch.prevGas;
        }
    }
    std::sort(res.selectedTxIds.begin(), res.selectedTxIds.end());
    return res;
}

BlockSelectionResult solveBlockGreedy(const std::vector<Tx>& txs, uint64_t gasLimit) {
    // Greedy fee-per-gas picking available nonces
    std::map<int, std::vector<Tx>> byAccount;
    for (const auto& tx : txs) {
        byAccount[tx.account].push_back(tx);
    }
    for (auto& pair : byAccount) {
        std::sort(pair.second.begin(), pair.second.end(), [](const Tx& a, const Tx& b) {
            return a.nonce < b.nonce;
        });
    }

    std::map<int, size_t> nextIndex;
    for (const auto& pair : byAccount) {
        nextIndex[pair.first] = 0;
    }

    BlockSelectionResult res;
    uint64_t remainingGas = gasLimit;

    while (true) {
        double bestRatio = -1.0;
        int bestAcc = -1;

        for (const auto& pair : byAccount) {
            int acc = pair.first;
            size_t idx = nextIndex[acc];
            if (idx < pair.second.size()) {
                const auto& tx = pair.second[idx];
                if (tx.gas <= remainingGas) {
                    double ratio = (double)tx.fee / (double)tx.gas;
                    if (ratio > bestRatio) {
                        bestRatio = ratio;
                        bestAcc = acc;
                    }
                }
            }
        }

        if (bestAcc < 0) break;

        const auto& chosenTx = byAccount[bestAcc][nextIndex[bestAcc]];
        res.selectedTxIds.push_back(chosenTx.id);
        res.totalFee += chosenTx.fee;
        res.totalGas += chosenTx.gas;
        remainingGas -= chosenTx.gas;
        nextIndex[bestAcc]++;
    }

    std::sort(res.selectedTxIds.begin(), res.selectedTxIds.end());
    return res;
}

}  // namespace merkle
