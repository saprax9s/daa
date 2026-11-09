# BASELINES.md: Exact Motif-Free Codon Design

Project: Exact Motif-Free Codon Design Using Aho-Corasick and Dynamic Programming (DAA course case study, single author).

This file lists every baseline and oracle: what it is, why it is included, its parameters, and its fairness notes. Hypotheses are expectations to be tested, not results. Results go in `docs/REPORT.md`, generated from `results/*.csv`.

## 0. Notation and fairness rules

- n = protein length. P = set of forbidden motifs. |Q| = number of automaton states. L = length of the longest expanded motif.
- [lo, hi] = global GC-count range. k, [Lw, Uw] = window size and window GC range.
- Score S = sum of per-codon scores plus adjacent codon-pair scores (integers, so comparisons are exact).

Rules that apply to every method:

1. Same instance files (`instances/*.json`), generated once by a fixed-seed generator and read by every implementation.
2. Same scores, same constraints (forbidden motifs on both strands, GC range, window GC), same codon table.
3. Every output is checked by one independent validity checker (motif search by naive substring match, GC counts, protein translation). A method's "validity rate" is judged by this checker only, never by the method itself.
4. Report failure types separately: "no valid sequence found" versus "valid but suboptimal".
5. Same timing method for C++ methods: release build, 5 warm-up runs, at least 30 repetitions, median with min and IQR, instance loading outside the timed region, automaton build timed separately.
6. Python baselines are never timed against C++. Report them for correctness and quality only.

## 1. Methods under comparison

| ID | Method | Exact? | Purpose |
|---|---|---|---|
| D | The DP with Aho-Corasick state (proposed method) | Yes | Subject of the study |
| M1 | Brute force / backtracking | Yes | Ground truth on small instances |
| M2 | Suffix-string oracle (Python) | Yes | Independent exact method for long proteins |
| M3 | Greedy | No | Shows when per-residue choices fail |
| M4 | Sample-and-reject | No | Shows the cost of generate-and-test |
| M5 | DNA Chisel (optional) | No (local search) | Widely used existing tool |
| M6 | 0-1 integer linear program (optional) | Yes | Independent exact comparison on moderate sizes |

### D: DP with automaton state (the method itself)

- **State:** (last codon c, automaton state q, GC count g). Window mode adds the GC counts of the last k - 1 codons.
- **Configurations to run as ablations:** D0 = no pruning; D1 = GC-bound pruning; D2 = D1 + dominance pruning (only valid when lo <= 0). Every configuration must return the same optimum score.
- **Complexity bound to compare against:** at most 6 * |Q| * (3i + 1) states at layer i, each with at most 6 transitions, so O(n^2 * |Q|) up to a constant of about 36 to 108.

### M1: Brute force / backtracking

- **What:** enumerate every codon sequence, test validity with the independent checker, keep the best score.
- **Why:** unquestionable ground truth.
- **Limits:** full genetic code only for n <= 5 (at most 6^n sequences). For n from 6 to 14 use the reduced code (exactly 2 codons per amino acid, 1 where only 1 exists, fixed deterministic choice).
- **Fairness:** same score function and constraints as D. Must report infeasibility when no valid sequence exists.

### M2: Suffix-string oracle (Python)

- **What:** an exact DP whose state is the last L - 1 nucleotides plus the GC count. No automaton.
- **Why:** an independent exact method for proteins too long for brute force. Agreement with D means the automaton logic is right.
- **Rules:** written without reading the C++ code. State count is bounded by 4^(L-1) times the GC range, so keep L modest (for example motifs up to 7 or 8 letters).
- **Fairness:** compared for correctness only. Not timed against C++.

### M3: Greedy

- **G1 (pure):** for each residue in order, choose the highest-scoring codon that creates no forbidden motif at that step. Deterministic tie-breaking (lowest codon index). GC range checked at the end.
- **G2 (GC-steered):** as G1, but when several codons are valid, prefer the one that keeps the running GC on track towards the middle of [lo, hi].
- **Why:** shows that a choice that is locally best can be globally bad or infeasible.
- **Report:** fraction of instances with no valid sequence, fraction valid but suboptimal, and score gap to D.
- **Fairness:** both variants use the same codon scores as D. Neither sees future residues.

### M4: Sample-and-reject

- **What:** sample codon sequences at random with probability driven by the same scores, test validity with the independent checker, and keep the best valid one found within a budget.
- **Parameters:** trial budgets B in {10^3, 10^4, 10^5, 10^6}; at least 3 seeds per instance; the temperature of the score-based sampling.
- **Report:** success rate, best score versus D, number of wasted candidates, and a time-matched variant (same wall-clock as D).
- **Fairness:** same constraints and scores. Report both budget-matched and time-matched results. Wasted candidates grow quickly with constraint tightness (long proteins, many motifs, narrow GC range).

### M5 (optional): DNA Chisel

- **What:** the DNA Chisel Python package, a constraint-based local-search optimizer. Used with avoid-pattern constraints on each motif and its reverse complement, and global or window GC-content constraints.
- **Why:** it is the tool a practitioner would use today.
- **Differences to document:**
  - Its default objectives (for example codon usage) differ from our additive pair score. For comparison, run it with the constraints only and then evaluate its output with our score and the independent checker.
  - It is stochastic. Fix seeds and run at least 3 seeds.
  - Record the package version.
- **Fairness:** compare validity rate and score gap only. Not timed against C++.
- **Honest outcome rule:** if DNA Chisel matches D's score on every tested constraint set, report that. The contribution is then the analysis and the exactness guarantee.

### M6 (optional): 0-1 integer linear program

- **What:** one binary variable per residue and codon, one codon per residue, linear GC constraints, and linear constraints that forbid codon combinations that create a motif. Adjacent-pair scores need extra linearised variables. Follows the 0-1 formulation of Arbib, Rossi and Tessitore (2020). Check the exact formulation in the paper before implementing.
- **Solver:** OR-Tools CP-SAT or PuLP with CBC. Record the solver version and the time limit.
- **Why:** an independent exact method with a very different structure. Agreement confirms exactness. Runtime growth shows why DP is attractive.
- **Fairness:** report solved-to-optimality rate and timeouts. Python-based, so no timing comparison against C++.

## 2. Validation oracles and tests

| ID | Check | Pass condition |
|---|---|---|
| V1 | D versus M1 | Identical optimum score. Infeasibility reported exactly when M1 finds nothing. |
| V2 | Automaton versus naive substring search | Identical verdicts on at least 100000 random DNA strings. Boundary-spanning motifs detected. |
| V3 | D versus M2 on long proteins; D versus M1 on the reduced code | Identical optimum. |
| V4 | M3 and M4 versus D | A measurable fraction of failures or suboptimal results. If none, report as a negative result. |
| V5 | Measured scaling versus the bound | Slope consistent with O(n^2 * |Q|). Pruning never changes the optimum. |
| V6 | Window-GC D versus M1 | Identical optimum. State count matches the counting bound. |

Any mismatch in V1, V2, V3 or V6 is a bug. Fix it before reporting any benchmark.

## 3. Instances and motif sets

- **Small:** random proteins with n up to 5, full genetic code (for M1).
- **Medium:** random proteins with n from 6 to 14, reduced code (for M1).
- **Long:** random proteins with n from about 20 to several hundred (for M2, and M6 where feasible).
- **Real proteins:** only if downloaded from a public database with the accession recorded. Never typed from memory.
- **Motif sets:**
  - Random short motifs (for tests).
  - Common restriction sites, both strands: EcoRI GAATTC, BamHI GGATCC, HindIII AAGCTT, XhoI CTCGAG, NdeI CATATG, NcoI CCATGG, XbaI TCTAGA, SpeI ACTAGT, PstI CTGCAG, KpnI GGTACC, SacI GAGCTC, NotI GCGGCCGC, BsaI GGTCTC, BsmBI CGTCTC, SapI GCTCTTC. Check each against REBASE where network access allows, and mark each entry as verified or unverified.
  - One IUPAC-ambiguous stress motif that expands to many concrete strings.
- **GC ranges:** loose and tight settings, plus a deliberately infeasible setting to test infeasibility reporting.

## 4. Metrics

- Optimality gap (score difference from D), validity rate, failure type, wasted candidates (M4).
- Runtime (C++ only), number of DP states, peak memory.
- Pruning ablation: states and time for D0, D1, D2.

## 5. Threats to validity

- Exactness holds only for the additive score model. Neural scores are not covered.
- No claim is made about expression or folding in a cell. Only validity and optimality on the stated objective are established.
- Synthetic scores may favour or disfavour a baseline. Where a real codon-usage table is used, it must come from a cited public source.
- Nucleotide-aligned windows are only handled conservatively. Codon-aligned windows are exact.
- Heavily ambiguous motifs enlarge the automaton and the DP table.
- Restriction sites not checked against REBASE are labelled unverified.
- Build and run times depend on the machine. Record the CPU, OS, compiler, and flags in `results/environment.txt`.

## 6. Intentionally not included

- LinearDesign-style folding-energy optimisation. This case study uses additive scores only.
- Neural language-model decoders.
- Any baseline whose objective cannot be evaluated under the shared additive score.
