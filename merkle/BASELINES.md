# BASELINES: Depth-Bounded Alphabetic Merkle Trees

Project: Depth-Bounded Alphabetic Merkle Trees for Stateless Blockchain Verification (DAA course project, single author).

This file lists every baseline: what it is, why it is included, its parameters, and its fairness notes. Hypotheses below are expectations to be tested, not results. Results go in `docs/REPORT.md`, generated from `results/*.csv`.

## 0. Notation and fairness rules

- N = number of accounts. K = number of hot accounts. The interleaved layout has 2K + 1 top-level leaves (K hot keys and K + 1 cold-bucket super-leaves).
- m = accounts touched by one block. s = Zipf exponent. Dmax = depth cap. L = number of top levels cached by a client.

Rules that apply to every baseline:

1. Same keys, same hash function (self-implemented SHA-256 with the same domain tags), same seeds, same workloads.
2. Same proof accounting: count hashes, bytes, and direction bits separately. Leaf data is counted separately.
3. Same multiproof algorithm for every tree (deduplicated shared hashes). No tree is compared against "m separate proofs" unless that column is labelled as such.
4. Same timing method: release build, 5 warm-up runs, at least 30 repetitions, median with min and IQR, workload generation outside the timed region.
5. Cold accounts are always reported, not only hot ones. A cold account's proof length is its top-tree depth plus its bucket subtree depth.
6. Every baseline reports hot-account proof size, cold-account proof size, and traffic-weighted overall size.

## 1. Tree-shape baselines

| ID | Baseline | Time | Space | Preserves key order | Why it is included |
|---|---|---|---|---|---|
| B1 | Full-state balanced binary Merkle tree over all N accounts | O(N) | O(N) | Yes | **Main baseline.** What blockchains use today: depth about log2 N for every account. |
| B2 | Balanced tree over the same 2K + 1 interleaved leaves (unweighted) | O(K) | O(K) | Yes | Isolates the effect of weighting. Same layout, weights ignored. |
| B3 | Mehlhorn weight-balanced tree | O(K) in theory; simple versions are O(K log K), state which one is implemented | O(K) | Yes | Fast heuristic. Known bound: cost within about 2 of the entropy. |
| B4 | Garsia-Wachs optimal alphabetic tree | O(K log K) if the full algorithm is implemented; state if a simpler O(K^2) variant is used | O(K) | Yes | Optimal without a depth cap. |
| B5 | Knuth interval DP with monotone split points | O(K^2) | O(K^2) | Yes | Second optimal method. Its cost must equal B4 exactly. |
| B6 | Depth-capped alphabetic DP, DP[h][i][j] | O(K^2 * Dmax) | O(K^2) with rolled layers and compact split storage | Yes | The proposed method. Hard cap on any leaf depth. |
| B7 | Weight-clamped Garsia-Wachs (every weight raised to at least W_total / 2^c) | O(K log K) | O(K) | Yes | Cheap defence against adversarial depth. Tests whether the DP cap is needed. |

**Parameters**

- K in {64, 128, 256, 512, 1024}.
- N in {2^16, 2^20}, fully hashed. N = 10^7 is analytic only (depth estimate, no full hashing).
- Dmax from the feasibility bound ceil(log2(2K + 1)) up to that bound plus 4.
- Clamp exponent c swept over a small range.

**Fairness notes**

- Weighted average depth must use identical leaf weights for B2 to B7.
- B4 and B5 are the same optimum. A cost mismatch means a bug.
- B6 equals B4 when the cap does not bind. Otherwise B6 cost is at least B4 cost.
- B1 versus B2 shows what the interleaved layout costs: cold accounts pay for their bucket subtree.

## 2. Related-work baselines

### B8: Huffman-layout hot tier with balanced cold tier (order NOT preserved)

- **What:** hot accounts in a Huffman tree (greedy merging by frequency), cold accounts in a balanced tree, combined under one root. Modelled on the idea of the Huffman-Merkle Tree in "Authenticated Data Structures for Dynamic Workloads" (IACR ePrint 2026/1235).
- **Why:** it is the closest recent work. It shows what is gained or lost by dropping sorted order.
- **Assumptions (document and keep):**
  - A1: re-implemented from a description, not from the authors' code. Verify against the full paper before citing any comparison as a faithful reproduction.
  - A2: membership proof size is measured directly on the tree.
  - A3: non-inclusion of a hot-tier key needs an auxiliary sorted structure. Model it as an extra balanced sorted tree over the hot keys, and report membership-only and membership + non-inclusion costs separately.
- **Fairness:** same weights, same K, same hash domain tags.

### B9: Static top-L-level caching

- **What:** the client caches the top L levels of the tree, so those hashes are not sent.
- **Why:** a cheap alternative to reshaping the tree.
- **Parameters:** L in {0, 2, 4, 6}.
- **Fairness:** apply the same caching rule to every tree (cache every node with depth < L), including the proposed tree. Otherwise the comparison is biased.

### B10 (optional): 16-ary Merkle tree

- **What:** wide Merkle tree with 16 children per node.
- **Why:** wide trees are a direction for Ethereum state, and they undercut the claim that balanced binary trees are the only alternative.
- **Sanity check (analytic):** one membership proof carries 15 sibling hashes per level. For N = 2^20 there are 5 levels, so 75 hashes, versus 20 for a binary tree.
- **Fairness:** report hashes and bytes. A single wide node costs more to verify, so also report verification time.

### B11 (optional): Simplified Merkle Patricia Trie

- **What:** hex-nibble trie over hashed keys with branch and extension nodes, without RLP encoding.
- **Why:** the structure Ethereum uses today.
- **Assumption:** simplified, so absolute byte sizes are approximate. Use it for proof-length comparison only and state this clearly.

## 3. Multiproof baseline

- **What:** deduplicated multiproof for a block that touches m accounts. Shared upper-level hashes are sent once.
- **Why:** comparing against m separate proofs would overstate any gain. Most of the saving comes from deduplication, which any tree gets.
- **Reference:** Compact Merkle Multiproofs (arXiv:2002.07648).
- **Parameters:** m in {20, 50, 100, 200, 500}; s in {0.8, 1.0, 1.15, 1.3}; cold-touch fraction swept (for example 0%, 10%, 30%).
- **Fairness:** same algorithm for every tree. Report hash count and bytes with direction bits and leaf data counted separately.

## 4. Adversarial-analysis baselines

| ID | Baseline | Purpose |
|---|---|---|
| A1 | Balanced tree (B1 or B2) | Immune to weight manipulation. Reference line for worst-case depth. |
| A2 | Unconstrained optimal trees (B3, B4, B5) | Expected to be vulnerable to geometric or Fibonacci weight cascades. |
| A3 | Weight-clamped tree (B7) | Does clamping give enough protection? |
| A4 | Depth-capped DP (B6) | Hard guarantee: max depth at most Dmax. |

- **Attack models:** (a) direct geometric or Fibonacci weights on 64 leaves; (b) budget-limited attacker, with per-block touches limited by a gas-derived cap and per-epoch totals limited by epoch length; (c) griefing, where the attacker pushes honest hot accounts deeper.
- **Metrics:** maximum depth reached versus attacker budget, minimum budget to reach depth d, and average-depth damage to honest accounts.
- **Fairness note:** a tree with a very low average depth under attack can still have a very long worst-case proof. Always report both.

## 5. Pre-block pipeline baselines (secondary scope)

| ID | Baseline | Purpose |
|---|---|---|
| P1 | Brute-force directed cycle checker (small graphs only) | Ground truth for Tarjan SCC. Results must agree. |
| P2 | Undirected union-find grouping | Negative control. Expected to flag every connected group, including harmless ones. Report its false-positive rate against Tarjan on DAG-shaped transfers. |
| P3 | Greedy fee-per-gas packing | Baseline for the knapsack. Report the actual gain. No assumed percentage. |
| P4 | Greedy packing with a two-resource fee density (gas and blob space) | Stronger greedy baseline. |
| P5 | Knapsack ignoring nonce order | Shows how many blocks become invalid if nonce prefixes are ignored. |
| P6 | Exhaustive search (small instances) | Ground truth for the nonce-prefix, gas-bucketed 2D knapsack DP. |

- **Gas bucketing:** 100k-gas buckets. State the rounding rule, since it affects optimality.

## 6. Comparison logic (what each pair isolates)

| Compare | What it shows |
|---|---|
| B1 vs B2 | Cost of the interleaved layout (cold accounts pay for their bucket subtree) |
| B2 vs B3 / B4 / B5 | Benefit of weighting |
| B3 vs B4 | Gap between the fast heuristic and the optimum |
| B4 vs B5 | Correctness cross-check (costs must match) |
| B4 vs B6 | Price of the depth cap in average depth |
| B6 vs B7 | Hard cap versus cheap weight clamping |
| B4 vs B8 | Cost and benefit of preserving sorted order |
| B9 on every tree | Whether caching top levels closes the gap without reshaping |
| Multiproof across trees | Real saving after deduplication |

## 7. Threats to validity

- Workloads are synthetic (Zipf). Real Ethereum traces are used only if they can be obtained, and are reported separately.
- B8 and B11 are re-implemented from descriptions. State every assumption.
- Hash-count comparisons ignore the cost of encoding direction bits and shape. Both are reported separately.
- Cold-bucket maintenance under account insertion is measured, not analysed in closed form.
- Build-time results depend on the machine. Record the CPU, OS, compiler, and flags in `results/environment.txt`.

## 8. Intentionally not included

- Zero-knowledge proof systems and Verkle commitments with polynomial openings, which are a different cost model.
- Any baseline that cannot support non-inclusion proofs without a stated extra cost. B8 is included only with the A3 assumption above.
