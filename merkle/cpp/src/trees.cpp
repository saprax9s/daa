#include "trees.h"

#include <algorithm>
#include <cstdint>
#include <functional>
#include <limits>
#include <vector>

namespace merkle {

namespace {

constexpr uint64_t INF = std::numeric_limits<uint64_t>::max() / 4;

// Balanced binary tree over sorted leaves, recursive midpoint split.
void balancedDepths(int i, int j, uint64_t d, std::vector<uint64_t>& depths) {
    if (j - i == 1) {
        depths[i] = d;
        return;
    }
    int m = (i + j) / 2;
    balancedDepths(i, m, d + 1, depths);
    balancedDepths(m, j, d + 1, depths);
}

}  // namespace

TreeResult solveBalanced(const std::vector<uint64_t>& w) {
    TreeResult r;
    r.leafCount = (uint64_t)w.size();
    std::vector<uint64_t> depths(w.size(), 0);
    if (!w.empty()) balancedDepths(0, (int)w.size(), 0, depths);
    uint64_t cost = 0;
    uint64_t maxd = 0;
    for (size_t i = 0; i < w.size(); ++i) {
        cost += w[i] * depths[i];
        maxd = std::max(maxd, depths[i]);
    }
    r.cost = cost;
    r.maxDepth = maxd;
    return r;
}

TreeResult solveOptimal(const std::vector<uint64_t>& w) {
    TreeResult r;
    r.leafCount = (uint64_t)w.size();
    int n = (int)w.size();
    if (n == 0) return r;
    if (n == 1) {
        r.cost = 0;
        r.maxDepth = 0;
        return r;
    }
    // Garsia-Wachs Phase 1
    struct GWNode {
        uint64_t weight;
        int leafIdx;
        int left = -1;
        int right = -1;
    };
    std::vector<GWNode> tree;
    tree.reserve(2 * n);
    std::vector<int> active;
    active.reserve(n);
    for (int i = 0; i < n; ++i) {
        GWNode node;
        node.weight = w[i];
        node.leafIdx = i;
        tree.push_back(node);
        active.push_back(i);
    }
    while (active.size() > 1) {
        int k = -1;
        for (size_t i = 1; i < active.size(); ++i) {
            uint64_t rw = (i + 1 < active.size()) ? tree[active[i + 1]].weight : INF;
            if (tree[active[i - 1]].weight <= rw) {
                k = (int)i;
                break;
            }
        }
        if (k == -1) k = (int)active.size() - 1;
        int a = active[k - 1];
        int b = active[k];
        GWNode merged;
        merged.weight = tree[a].weight + tree[b].weight;
        merged.leafIdx = -1;
        merged.left = a;
        merged.right = b;
        int newId = (int)tree.size();
        tree.push_back(merged);
        active.erase(active.begin() + k);
        active.erase(active.begin() + (k - 1));
        int insertPos = 0;
        for (int j = k - 2; j >= 0; --j) {
            if (tree[active[j]].weight > merged.weight) {
                insertPos = j + 1;
                break;
            }
        }
        active.insert(active.begin() + insertPos, newId);
    }
    // Phase 2: compute depths
    std::vector<uint64_t> depths(n, 0);
    std::function<void(int, uint64_t)> getDepths = [&](int id, uint64_t d) {
        if (tree[id].leafIdx >= 0) {
            depths[tree[id].leafIdx] = d;
            return;
        }
        getDepths(tree[id].left, d + 1);
        getDepths(tree[id].right, d + 1);
    };
    getDepths(active[0], 0);

    uint64_t cost = 0;
    uint64_t maxd = 0;
    for (int i = 0; i < n; ++i) {
        cost += w[i] * depths[i];
        maxd = std::max(maxd, depths[i]);
    }
    r.cost = cost;
    r.maxDepth = maxd;
    return r;
}

TreeResult solveKnuth(const std::vector<uint64_t>& w) {
    TreeResult r;
    r.leafCount = (uint64_t)w.size();
    int n = (int)w.size();
    if (n == 0) return r;
    if (n == 1) {
        r.cost = 0;
        r.maxDepth = 0;
        return r;
    }
    std::vector<uint64_t> pref(n + 1, 0);
    for (int i = 0; i < n; ++i) pref[i + 1] = pref[i] + w[i];
    auto W = [&](int i, int j) { return pref[j + 1] - pref[i]; };

    std::vector<std::vector<uint64_t>> C(n, std::vector<uint64_t>(n, 0));
    std::vector<std::vector<int>> opt(n, std::vector<int>(n, -1));
    for (int i = 0; i < n; ++i) {
        C[i][i] = 0;
        opt[i][i] = i;
    }
    for (int len = 2; len <= n; ++len) {
        for (int i = 0; i + len - 1 < n; ++i) {
            int j = i + len - 1;
            uint64_t best = INF;
            int bestm = -1;
            int lo = (len == 2) ? i : opt[i][j - 1];
            int hi = (len == 2) ? i : opt[i + 1][j];
            for (int m = lo; m <= hi && m < j; ++m) {
                uint64_t c = C[i][m] + C[m + 1][j] + W(i, j);
                if (c < best) {
                    best = c;
                    bestm = m;
                }
            }
            C[i][j] = best;
            opt[i][j] = bestm;
        }
    }
    std::vector<uint64_t> depths(n, 0);
    std::function<void(int, int, uint64_t)> rec = [&](int i, int j, uint64_t d) {
        if (i == j) {
            depths[i] = d;
            return;
        }
        int m = opt[i][j];
        rec(i, m, d + 1);
        rec(m + 1, j, d + 1);
    };
    rec(0, n - 1, 0);
    uint64_t cost = 0;
    uint64_t maxd = 0;
    for (int i = 0; i < n; ++i) {
        cost += w[i] * depths[i];
        maxd = std::max(maxd, depths[i]);
    }
    r.cost = cost;
    r.maxDepth = maxd;
    return r;
}

TreeResult solveDepthCapped(const std::vector<uint64_t>& w, uint64_t Dmax) {
    TreeResult r;
    r.leafCount = (uint64_t)w.size();
    int n = (int)w.size();
    if (n == 0) return r;
    if (n == 1) {
        r.cost = 0;
        r.maxDepth = 0;
        return r;
    }
    if (Dmax == 0) {
        r.cost = INF;
        r.maxDepth = 1;
        return r;
    }
    uint64_t capPower = 1ULL;
    bool overflow = false;
    for (uint64_t i = 0; i < Dmax; ++i) {
        capPower *= 2ULL;
        if (capPower >= (uint64_t)n) {
            overflow = true;
            break;
        }
    }
    if (!overflow && capPower < (uint64_t)n) {
        r.cost = INF;
        r.maxDepth = Dmax + 1;
        return r;
    }

    std::vector<uint64_t> pref(n + 1, 0);
    for (int i = 0; i < n; ++i) pref[i + 1] = pref[i] + w[i];
    auto W = [&](int i, int j) { return pref[j + 1] - pref[i]; };

    std::vector<std::vector<std::vector<uint64_t>>> dp(
        Dmax + 1, std::vector<std::vector<uint64_t>>(n, std::vector<uint64_t>(n, INF)));
    std::vector<std::vector<std::vector<int16_t>>> split(
        Dmax + 1, std::vector<std::vector<int16_t>>(n, std::vector<int16_t>(n, -1)));

    for (int i = 0; i < n; ++i) dp[0][i][i] = 0;

    for (int h = 1; h <= (int)Dmax; ++h) {
        for (int i = 0; i < n; ++i) dp[h][i][i] = 0;
        uint64_t maxLeavesForH = (h >= 60) ? INF : (1ULL << h);
        uint64_t maxSubLeaves = (h - 1 >= 60) ? INF : (1ULL << (h - 1));
        for (int len = 2; len <= n; ++len) {
            if ((uint64_t)len > maxLeavesForH) break;
            for (int i = 0; i + len - 1 < n; ++i) {
                int j = i + len - 1;
                uint64_t best = INF;
                int16_t bestm = -1;
                int lo = (maxSubLeaves >= (uint64_t)len) ? i : (int)(j - maxSubLeaves);
                int hi = (maxSubLeaves >= (uint64_t)len) ? j - 1 : (int)(i + maxSubLeaves - 1);
                lo = std::max(i, lo);
                hi = std::min(j - 1, hi);
                for (int m = lo; m <= hi; ++m) {
                    if (dp[h - 1][i][m] >= INF || dp[h - 1][m + 1][j] >= INF) continue;
                    uint64_t c = dp[h - 1][i][m] + dp[h - 1][m + 1][j] + W(i, j);
                    if (c < best) {
                        best = c;
                        bestm = (int16_t)m;
                    }
                }
                dp[h][i][j] = best;
                split[h][i][j] = bestm;
            }
        }
    }

    if (dp[Dmax][0][n - 1] >= INF) {
        r.cost = INF;
        r.maxDepth = Dmax + 1;
        return r;
    }

    std::vector<uint64_t> depths(n, 0);
    std::function<void(int, int, int, uint64_t)> rec = [&](int i, int j, int h, uint64_t d) {
        if (i == j) {
            depths[i] = d;
            return;
        }
        int m = split[h][i][j];
        if (m < i || m >= j) m = (i + j) / 2;
        rec(i, m, h - 1, d + 1);
        rec(m + 1, j, h - 1, d + 1);
    };
    rec(0, n - 1, (int)Dmax, 0);

    uint64_t cost = 0;
    uint64_t maxd = 0;
    for (int i = 0; i < n; ++i) {
        cost += w[i] * depths[i];
        maxd = std::max(maxd, depths[i]);
    }
    r.cost = cost;
    r.maxDepth = maxd;
    return r;
}

TreeResult solveMehlhorn(const std::vector<uint64_t>& w) {
    TreeResult r;
    r.leafCount = (uint64_t)w.size();
    int n = (int)w.size();
    if (n == 0) return r;
    if (n == 1) {
        r.cost = 0;
        r.maxDepth = 0;
        return r;
    }
    std::vector<uint64_t> pref(n + 1, 0);
    for (int i = 0; i < n; ++i) pref[i + 1] = pref[i] + w[i];
    auto W = [&](int i, int j) { return pref[j + 1] - pref[i]; };

    std::vector<uint64_t> depths(n, 0);
    std::function<void(int, int, uint64_t)> rec = [&](int i, int j, uint64_t d) {
        if (i == j) {
            depths[i] = d;
            return;
        }
        uint64_t total = W(i, j);
        int bestm = i;
        uint64_t bestDiff = INF;
        for (int m = i; m < j; ++m) {
            uint64_t leftW = W(i, m);
            uint64_t diff = (2 * leftW >= total) ? (2 * leftW - total) : (total - 2 * leftW);
            if (diff < bestDiff) {
                bestDiff = diff;
                bestm = m;
            }
        }
        rec(i, bestm, d + 1);
        rec(bestm + 1, j, d + 1);
    };
    rec(0, n - 1, 0);
    uint64_t cost = 0;
    uint64_t maxd = 0;
    for (int i = 0; i < n; ++i) {
        cost += w[i] * depths[i];
        maxd = std::max(maxd, depths[i]);
    }
    r.cost = cost;
    r.maxDepth = maxd;
    return r;
}

TreeResult solveClamped(const std::vector<uint64_t>& w, int c) {
    TreeResult r;
    r.leafCount = (uint64_t)w.size();
    if (w.empty()) return r;
    uint64_t total = 0;
    for (uint64_t x : w) total += x;
    uint64_t denom = 1ULL << c;
    uint64_t threshold = (total + denom - 1) / denom;
    std::vector<uint64_t> cw = w;
    for (auto& x : cw) {
        if (x < threshold) x = threshold;
    }
    return solveOptimal(cw);
}

TreeResult solveHuffman(const std::vector<uint64_t>& hotWeights, const std::vector<uint64_t>& coldBucketWeights) {
    TreeResult r;
    std::vector<uint64_t> hw = hotWeights;
    std::vector<std::vector<int>> groups;
    std::vector<uint64_t> wts;
    for (size_t i = 0; i < hw.size(); ++i) {
        std::vector<int> g(1, (int)i);
        groups.push_back(g);
        wts.push_back(hw[i]);
    }
    std::vector<int> depthOf(hw.size(), 0);
    while (groups.size() > 1) {
        int bi = 0;
        for (size_t i = 1; i < groups.size(); ++i) {
            if (wts[i] < wts[bi]) bi = (int)i;
        }
        int bj = -1;
        for (size_t i = 0; i < groups.size(); ++i) {
            if ((int)i == bi) continue;
            if (bj < 0 || wts[i] < wts[bj]) bj = (int)i;
        }
        uint64_t nw = wts[bi] + wts[bj];
        std::vector<int> g = groups[bi];
        for (int x : groups[bj]) {
            depthOf[x] += 1;
            g.push_back(x);
        }
        for (int x : groups[bi]) depthOf[x] += 1;
        groups[bi] = g;
        wts[bi] = nw;
        groups.erase(groups.begin() + bj);
        wts.erase(wts.begin() + bj);
    }
    uint64_t huffCost = 0;
    uint64_t huffMax = 0;
    for (size_t i = 0; i < hw.size(); ++i) {
        huffCost += hw[i] * (uint64_t)depthOf[i];
        huffMax = std::max(huffMax, (uint64_t)depthOf[i]);
    }
    TreeResult cr = solveBalanced(coldBucketWeights);
    r.cost = huffCost + cr.cost;
    r.maxDepth = std::max(huffMax, cr.maxDepth);
    r.leafCount = (uint64_t)(hw.size() + coldBucketWeights.size());
    return r;
}

uint64_t cachedCost(const std::vector<uint64_t>& w, uint64_t L, const std::vector<uint64_t>& depths) {
    uint64_t cost = 0;
    for (size_t i = 0; i < w.size(); ++i) {
        uint64_t d = depths[i];
        if (d > L) cost += w[i] * (d - L);
    }
    return cost;
}

MerkleTree buildOptimalTree(const std::vector<uint64_t>& w) {
    MerkleTree t;
    int n = (int)w.size();
    if (n == 0) return t;
    if (n == 1) {
        TreeNode leaf;
        leaf.isLeaf = true;
        leaf.leafIdx = 0;
        leaf.left = -1;
        leaf.right = -1;
        leaf.parent = -1;
        t.nodes.push_back(leaf);
        t.root = 0;
        t.depths = {0};
        t.hashes.resize(1);
        t.hashes[0] = hashAccount(0, w[0]);
        t.rootHash = t.hashes[0];
        return t;
    }
    std::vector<uint64_t> pref(n + 1, 0);
    for (int i = 0; i < n; ++i) pref[i + 1] = pref[i] + w[i];
    auto W = [&](int i, int j) { return pref[j + 1] - pref[i]; };
    std::vector<std::vector<uint64_t>> C(n, std::vector<uint64_t>(n, 0));
    std::vector<std::vector<int>> opt(n, std::vector<int>(n, -1));
    for (int i = 0; i < n; ++i) {
        C[i][i] = 0;
        opt[i][i] = i;
    }
    for (int len = 2; len <= n; ++len) {
        for (int i = 0; i + len - 1 < n; ++i) {
            int j = i + len - 1;
            uint64_t best = INF;
            int bestm = -1;
            int lo = (len == 2) ? i : opt[i][j - 1];
            int hi = (len == 2) ? i : opt[i + 1][j];
            for (int m = lo; m <= hi && m < j; ++m) {
                uint64_t c = C[i][m] + C[m + 1][j] + W(i, j);
                if (c < best) {
                    best = c;
                    bestm = m;
                }
            }
            C[i][j] = best;
            opt[i][j] = bestm;
        }
    }
    t.nodes.resize(n);
    for (int i = 0; i < n; ++i) {
        t.nodes[i].isLeaf = true;
        t.nodes[i].leafIdx = i;
        t.nodes[i].left = -1;
        t.nodes[i].right = -1;
        t.nodes[i].parent = -1;
    }
    std::function<int(int, int)> build = [&](int i, int j) -> int {
        if (i == j) return i;
        int m = opt[i][j];
        int leftChild = build(i, m);
        int rightChild = build(m + 1, j);
        TreeNode internalNode;
        internalNode.isLeaf = false;
        internalNode.leafIdx = -1;
        internalNode.left = leftChild;
        internalNode.right = rightChild;
        internalNode.parent = -1;
        int myId = (int)t.nodes.size();
        t.nodes.push_back(internalNode);
        t.nodes[leftChild].parent = myId;
        t.nodes[rightChild].parent = myId;
        return myId;
    };
    t.root = build(0, n - 1);
    t.depths.assign(n, 0);
    std::function<void(int, uint64_t)> computeDepths = [&](int id, uint64_t d) {
        if (t.nodes[id].isLeaf) {
            t.depths[t.nodes[id].leafIdx] = d;
            return;
        }
        computeDepths(t.nodes[id].left, d + 1);
        computeDepths(t.nodes[id].right, d + 1);
    };
    computeDepths(t.root, 0);
    t.hashes.resize(t.nodes.size());
    std::function<void(int)> fillHashes = [&](int id) {
        if (t.nodes[id].isLeaf) {
            t.hashes[id] = hashAccount((uint64_t)t.nodes[id].leafIdx, w[t.nodes[id].leafIdx]);
            return;
        }
        fillHashes(t.nodes[id].left);
        fillHashes(t.nodes[id].right);
        t.hashes[id] = hashInternal(t.hashes[t.nodes[id].left], t.hashes[t.nodes[id].right]);
    };
    fillHashes(t.root);
    t.rootHash = t.hashes[t.root];
    return t;
}

}  // namespace merkle
