ROLE
You are a senior algorithms/systems engineer. Implement, test, benchmark and document the project below. Be rigorous. Never fabricate numbers: every number in any document must come from a result file you generated. If you cannot execute code in this environment, write everything, leave clearly marked TBD placeholders for numbers, and give me exact commands to run.

PROJECT
"Depth-Bounded Alphabetic Merkle Trees for Stateless Blockchain Verification" (DAA course project, single author). Blockchain account traffic is Zipf-skewed. I compare Merkle tree shapes that put hot accounts near the root while keeping leaves sorted by address (needed for non-inclusion proofs). I do NOT claim new algorithms; this is an empirical study of known ones plus a depth cap and an adversarial analysis.

Base work (for context and baselines): classical optimal alphabetic trees (Hu-Tucker 1971, Garsia-Wachs 1977, Knuth 1971 DP, Mehlhorn 1975/77, height-restricted DP Wessner 1976 / Larmore 1987); recent: Huffman-Merkle Tree "Authenticated Data Structures for Dynamic Workloads" (IACR ePrint 2026/1235), traffic-aware Merkle trees (ePrint 2024/1451), QMDB (arXiv 2501.05262), AlDBaran (arXiv 2508.10493); Compact Merkle Multiproofs (arXiv 2002.07648); RFC 6962 hashing.

LANGUAGES
- C++20 is PRIMARY (all experiments, all timing). No external libraries (implement SHA-256 yourself). Must build with g++ (MinGW/Linux) and MSVC; provide a CMakeLists.txt plus simple build scripts (build.bat, build.sh).
- Rust is SECONDARY: port the core (SHA-256, scoring, tree constructions, proofs, a reduced benchmark) with no external crates except optional criterion. Rust output must use the same CSV schema as C++. If cargo is not installed, skip Rust gracefully and tell me how to install it (rustup) in the README.
- Python is allowed ONLY for plotting and report-table generation (matplotlib, pandas) and a tiny brute-force oracle for tiny inputs. Never time Python implementations.
- Identical PRNG (splitmix64) and identical workload generators in C++ and Rust so both produce identical roots and costs for the same seed. Add a cross-language equality test.

PART 1: REPOSITORY LAYOUT
repo/
  README.md              (exact run instructions, Windows PowerShell + Linux)
  BASELINES.md           (every baseline: what it is, why included, parameters, fairness notes)
  docs/REPORT.md         (living report, see Part 11)
  cpp/  (src/, tests/, bench/, CMakeLists.txt, build.bat, build.sh)
  rust/ (Cargo project)
  python/ (plots.py, make_tables.py)
  results/ (raw CSVs)   figures/
Single CSV schema: language,experiment,algo,K,N,m,zipf_s,dmax,seed,rep,metric,value,unit.

PART 2: PRIMITIVES
- SHA-256 from scratch. Test vectors: SHA256("abc") = ba7816bf 8f01cfea 414140de 5dae2223 b00361a3 96177a9c b410ff61 f20015ad; SHA256("") = e3b0c442 98fc1c14 9afbf4c8 996fb924 27ae41e4 649b934c a495991b 7852b855.
- Domain-separated hashing with a 1-byte tag: 0x00 account leaf, 0x01 internal node, 0x02 cold-bucket range root, 0x03 top-level state+weight commitment.
- Deterministic scoring, uint64 only, no floating point: W_new(a) = (3*touches(a)*1000000 + 7*W_old(a)) / 10 (integer division), updated once per epoch from touches only. Hot set = top K scores, ties broken by lowest account index. Commit a digest of the weight vector in the top-level commitment. Test: same inputs give the same root across runs, and across C++/Rust.
- Keys: random 64-bit addresses from the PRNG, sorted. Hot keys are scattered randomly in key space.

PART 3: TREE CONSTRUCTIONS (all take sorted leaf weights; all output a tree shape plus cost and max depth)
Interleaved layout: K sorted hot keys split the other N-K accounts into K+1 cold buckets. Each bucket is a balanced subtree (range-root tag 0x02) used as ONE super-leaf between its hot neighbours, giving a single sorted tree with 2K+1 top-level leaves. Bucket weight = sum of cold touches in its interval; every leaf has minimum weight 1.
Implement:
 1. Balanced (same 2K+1 leaves, unweighted)
 2. Mehlhorn weight-balanced (report the complexity you actually achieve)
 3. Garsia-Wachs optimal alphabetic. Prefer a true O(K log K) implementation; if you implement a simpler O(K^2) variant, say so explicitly in docs and in the measured slope.
 4. Knuth DP, O(K^2) via monotone split-point optimisation
 5. Depth-capped alphabetic DP, DP[h][i][j], O(K^2 * Dmax), requires Dmax >= ceil(log2(2K+1)). Control memory: roll the height layers, store split choices compactly (uint16), report peak memory.
 6. Weight-clamped Garsia-Wachs: every weight raised to at least W_total / 2^c (a cheap defence baseline against the cap)
Also report which of these produce identical optimal cost (3 vs 4 must match exactly).

PART 4: PROOFS
- Membership proof with direction bits (the shape varies per leaf).
- Non-inclusion proof: two adjacent leaves with a verifier that checks adjacency from the paths. If the key falls in a cold bucket, the proof extends into the bucket subtree.
- Batch multiproof with shared-hash deduplication.
- An independent verifier code path. Report proof size in hashes AND bytes, with direction bits and leaf data counted separately.

PART 5: WORKLOADS AND ADVERSARY
- Zipf account access generator, exponent s in {0.8, 1.0, 1.15, 1.3}, deterministic by seed. Blocks of T transactions, each touching 2 accounts; epochs of E blocks (default E=100).
- N in {2^16, 2^20} fully hashed; N=10^7 optional with analytic depth only (no full hashing).
- Adversary experiment: (a) direct geometric/Fibonacci weight cascades on 64 leaves; (b) budget-limited attacker: per-block touch cap derived from a gas limit, per-epoch total bounded; report the maximum depth each construction reaches vs attacker budget, and the minimum budget needed to reach depth d. Also a griefing test (attacker pushes honest hot accounts deeper).

PART 6: BASELINES (document all in BASELINES.md)
 B1 full-state balanced binary Merkle tree over all N accounts (MAIN baseline, depth ~log2 N)
 B2 balanced version of the same interleaved layout (isolates the effect of weighting)
 B3 Mehlhorn, B4 Garsia-Wachs, B5 Knuth, B6 depth-capped DP, B7 weight-clamped GW
 B8 Huffman-layout hot tier + balanced cold tier, in the style of the 2026 Huffman-Merkle Tree paper (ORDER NOT PRESERVED): measure membership proof size, and the extra cost of making non-inclusion proofs work. State in BASELINES.md that it is re-implemented from a description, with your assumptions.
 B9 static caching of the top L levels of the balanced tree
 B10 (optional) 16-ary Merkle tree, B11 (optional) simplified Merkle Patricia Trie
 Multiproof baseline: deduplicated balanced multiproof (not m separate proofs).

PART 7: EXPERIMENTS AND MEASUREMENT METHOD
E1 build time and operation counts (DP transitions, hash calls) vs K in {64,128,256,512,1024}
E2 proof size: hot accounts, cold accounts, traffic-weighted overall, against B1 on the full N-account tree
E3 adversarial depth (Part 5)
E4 block witness size vs m in {20,50,100,200,500}, vs s, vs cold-touch fraction
E5 epoch rebuild time and cold-bucket split/merge cost when accounts enter or leave the hot set
E6 Dmax sweep: average depth vs worst-case cap, with the feasibility bound
E7 non-inclusion proof size (all trees that support it)
E8 verification time per proof
E9 (secondary) pre-block pipeline: Tarjan SCC wash-trade detection O(V+E) on synthetic transfer graphs (with a brute-force cycle checker for small graphs), and nonce-prefix gas-bucketed 2D knapsack DP (compare with greedy fee-per-gas; report actual gain, not an assumed 30%).
Timing method: steady_clock (C++) / Instant (Rust); release builds (-O2 / --release); 5 warm-up runs, then at least 30 repetitions; report median, min, IQR; use sinks so the compiler cannot delete work; exclude I/O and workload generation from timed regions; fixed seeds; time hashing separately from tree construction; record CPU model, OS, compiler version and flags into results/environment.txt.

PART 8: TESTS (must pass before benchmarks)
SHA vectors; Garsia-Wachs cost == Knuth cost == brute-force optimum for K <= 10 (random weights); capped DP cost >= unconstrained cost, equal when the cap is not binding, and max depth <= Dmax always; Mehlhorn cost <= H + 2 (entropy bound); every valid membership/non-inclusion/multiproof verifies, tampered ones are rejected; leaf order preserved; determinism; C++/Rust equality; Tarjan vs brute force.

PART 9: COMPLEXITY ANALYSIS
For each algorithm give the theoretical time and space complexity, then fit the measured slope on a log-log plot (doubling ratio) and report theory vs measured in one table. Explain any mismatch.

PART 10: RUST PORT
Port after the C++ tests pass. Same CSV output, reduced experiments (E1, E2, E4). Compare C++ vs Rust runtimes on the same machine, and say clearly that they are not a language shootout.

PART 11: DOCUMENTATION, WRITTEN IN PARALLEL
Maintain docs/REPORT.md from the start and update it after every part (not at the end). Sections: 1 Problem statement; 2 Base work and what is reused; 3 Approach; 4 Algorithms and complexity (Part 9 table); 5 Experimental setup (machine, flags, seeds, workloads, baselines); 6 Results (tables and figures generated by python/make_tables.py from results/*.csv, never typed by hand); 7 Analysis (when each tree wins and loses, including cold-account cost, attack cost, depth-cap trade-off); 8 Limitations; 9 Reproduction. Report negative results honestly (for example if cold accounts have longer proofs than balanced, or the cap costs too much). No novelty claims.

PART 12: FINAL HAND-OFF (your last message must contain)
1) exact build and run commands for Windows (PowerShell, VS Code) and Linux, including the one-command run_all and expected runtime; 2) the list of baselines; 3) where each result file and figure lands; 4) how execution time is measured; 5) which tests passed; 6) anything unverified, simplified, or assumed.

WORK ORDER: Part 1 and 2, then 3 and 4, then 8 (tests), then 5 to 7, then 9 and 11, then 10 and E9. Stop after each stage, summarise what was done, and update the report.