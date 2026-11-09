# Exact Motif-Free Codon Design Using Aho-Corasick and Dynamic Programming

Design and Analysis of Algorithms (DAA) case study: Exact, polynomial-time reverse-translation of protein sequences into synonymous DNA sequences avoiding all forbidden restriction sites and motifs, satisfying global and sliding-window GC constraints, and maximizing an additive codon/pair-adaptation score.

---

## Repository Layout

```
codon/
├── README.md               # Quickstart and execution guide (Windows + Linux)
├── BASELINES.md            # Baseline descriptions (M1–M6) and fairness rules
├── docs/
│   └── REPORT.md           # Comprehensive living empirical report with real numbers
├── cpp/                    # Primary implementation (C++20, no external dependencies)
│   ├── CMakeLists.txt      # Cross-platform CMake configuration
│   ├── build.bat           # Windows PowerShell build script (direct g++ fallback)
│   ├── build.sh            # Linux / macOS build script
│   ├── src/                # genetic_code, automaton, dp, motifs, prng, scores, baselines
│   ├── tests/              # Selftest suite validating V1 through V6
│   └── bench/              # Benchmark harness for experiments E1-E8
├── python/                 # Independent verification oracle & visualization
│   ├── oracle.py           # Independent Python suffix-string DP oracle
│   ├── plots.py            # Generates figures/fig1..fig7.png
│   └── make_tables.py      # Formats Markdown tables from results/codon_results.csv
├── rust/                   # Optional secondary Rust port
│   ├── Cargo.toml
│   └── src/main.rs
├── instances/              # Canonical test and benchmark instances (JSON)
│   ├── real_gfp.json       # Aequorea victoria GFP (UniProt P42212, 238 aa)
│   ├── real_insulin.json   # Mature human insulin (UniProt P01308, 51 aa)
│   ├── small_full_n5.json
│   ├── small_reduced_n8.json
│   └── small_reduced_n12.json
├── results/
│   ├── environment.txt     # CPU model, OS, compiler version, optimization flags
│   └── codon_results.csv   # Raw empirical benchmark measurements (1283 rows)
└── figures/                # Publication-quality 300 DPI evaluation plots
    ├── fig1_optimality_gap.png
    ├── fig2_runtime_vs_n.png
    ├── fig3_runtime_vs_Q.png
    ├── fig4_pruning_ablation.png
    ├── fig5_window_gc.png
    ├── fig6_realistic_proteins.png
    └── fig7_sample_wasted.png
```

---

## Quickstart & Reproduction

### Prerequisites
- **C++ Compiler**: GCC 11+ / Clang 13+ / MSVC 2022 supporting `-std=c++20`.
- **Python**: Python 3.10+ with `pandas`, `numpy`, and `matplotlib`.
- **Rust (Optional)**: `cargo` 1.70+. See [Installing Rust](#installing-rust) below.

### One-Command Run (Windows PowerShell)

```powershell
# 1. Compile C++ code
cd M:\DAA\codon\cpp
.\build.bat

# 2. Run self-tests (Validations V1 through V6)
.\build\codon_tool.exe --selftest

# 3. Generate instances
.\build\codon_tool.exe --generate-instances M:\DAA\codon\instances

# 4. Run all experiments E1 to E8 (writes to results/codon_results.csv)
.\build\codon_tool.exe --bench M:\DAA\codon\results\codon_results.csv

# 5. Generate all figures and summary tables
cd ..
& "C:\Users\Sapra\AppData\Local\Programs\Python\Python313\python.exe" python/plots.py
& "C:\Users\Sapra\AppData\Local\Programs\Python\Python313\python.exe" python/make_tables.py
```

*Expected Execution Time*: **~12 seconds total**.

### One-Command Run (Linux / macOS)

```bash
cd codon/cpp
chmod +x build.sh
./build.sh
./build/codon_tool --selftest
./build/codon_tool --generate-instances ../instances
./build/codon_tool --bench ../results/codon_results.csv
cd ..
python3 python/plots.py
python3 python/make_tables.py
```

---

## Installing Rust

The primary implementation is in standard C++20 with no external libraries. A clean standalone Rust port is provided in `rust/`. If `cargo` is missing, install it via `rustup`:

- **Windows**:
  ```powershell
  winget install --id Rustlang.Rustup
  # or from https://rustup.rs
  ```
- **Linux / macOS**:
  ```bash
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
  source $HOME/.cargo/env
  ```

Once installed:
```bash
cd codon/rust
cargo run --release -- --test
cargo run --release -- --bench ../results/codon_rust_results.csv
```

---

## Key Findings

1. **Complexity Verification**: The measured DP runtime scales as $O(n^{2.06} \cdot |Q|)$, matching the theoretical $O(n^2 \cdot |Q|)$ bound within 3% error.
2. **Pruning Efficiency**: Dominance pruning (D2) shrinks the explored state space by **65.4%** and delivers a **$3.37\times$ speedup** with zero loss of optimality.
3. **Biological Gene Design**:
   - Mature Insulin (51 aa): Designed in **2.28 milliseconds** (6,389 states).
   - Full-length GFP (238 aa, 714 bp): Designed in **49.8 milliseconds** (123,514 states), avoiding all 15 Type II and Type IIS restriction sites across both strands.
4. **Heuristic Failure**: Greedy and random sampling fail to find any feasible sequence on 10%–15% of test instances and suffer a 60%–66% optimality penalty when feasible.
