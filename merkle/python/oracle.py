"""
Independent brute-force and DP oracle for Alphabetic Merkle Trees.
Written independently as ground truth for verifying C++ implementations.
"""

import sys
import math

def brute_optimal_alphabetic(weights):
    """O(N^3) standard alphabetic tree DP (Knuth / Gilbert-Moore)."""
    n = len(weights)
    if n == 0:
        return 0, 0
    if n == 1:
        return 0, 0
    pref = [0] * (n + 1)
    for i in range(n):
        pref[i + 1] = pref[i] + weights[i]
    W = lambda i, j: pref[j + 1] - pref[i]

    INF = float('inf')
    dp = [[0] * n for _ in range(n)]
    split = [[-1] * n for _ in range(n)]

    for length in range(2, n + 1):
        for i in range(n - length + 1):
            j = i + length - 1
            best = INF
            best_m = -1
            for m in range(i, j):
                c = dp[i][m] + dp[m + 1][j] + W(i, j)
                if c < best:
                    best = c
                    best_m = m
            dp[i][j] = best
            split[i][j] = best_m

    depths = [0] * n
    def get_depths(i, j, d):
        if i == j:
            depths[i] = d
            return
        m = split[i][j]
        get_depths(i, m, d + 1)
        get_depths(m + 1, j, d + 1)

    get_depths(0, n - 1, 0)
    cost = sum(weights[i] * depths[i] for i in range(n))
    return cost, max(depths)

def main():
    test_weights = [
        [1, 1, 1, 1],
        [2, 3, 1, 4, 5, 9, 2, 6],
        [10, 1, 1, 10],
        [5, 2, 8, 3, 1]
    ]
    print("Running Python independent oracle checks:")
    for w in test_weights:
        cost, max_d = brute_optimal_alphabetic(w)
        print(f"Weights: {w} -> Optimal Cost: {cost}, Max Depth: {max_d}")

if __name__ == "__main__":
    main()

