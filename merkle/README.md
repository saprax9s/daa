# Depth-Bounded Alphabetic Merkle Trees for Stateless Blockchain Verification

Implementation and empirical evaluation of **Depth-Bounded Alphabetic Merkle Trees** for stateless blockchain clients (Design and Analysis of Algorithms course case study, single author).

## Overview

Stateless blockchain verifiers require proof of account inclusion and non-inclusion to validate transactions without storing the full global state. Account accesses follow a heavy-tailed Zipf distribution ($s \ge 1.0$), where hot accounts dominate network traffic. This repository implements and benchmarks alphabetic Merkle tree structures that:
1. Preserve strict sorted key order (required for non-inclusion proofs).
2. Minimize proof length for high-frequency accounts using an interleaved super-leaf layout.
3. Impose a hard depth ceiling $D_{\max}$ to prevent adversarial proof-inflation attacks.

---

## Repository Structure

```
merkle/
├── README.md               # Quickstart and execution guide (Windows + Linux)
├── BASELINES.md            # Detailed specification of all 11 baselines and fairness rules
├── docs/
│   └── REPORT.md           # Living empirical report populated with genuine benchmark data
├── cpp/                    # Primary implementation (C++20, no external libraries)
│   ├── CMakeLists.txt      # Cross-platform CMake configuration
│   ├── build.bat           # Windows PowerShell build script (direct g++ fallback)
│   ├── build.sh            # Linux / macOS shell build script
│   ├── src/                # SHA-256, hashing, weights, trees, proofs, workloads, pipeline
│   ├── tests/              # Comprehensive selftest suite (17 checks)
│   └── bench/              # Benchmark harness for experiments E1-E9
├── python/                 # Independent verification oracle & visualization
│   ├── oracle.py           # Independent Python DP oracle for ground truth
│   ├── plots.py            # Generates figures/fig1..fig6.png
│   └── make_tables.py      # Extracts formatted Markdown tables from raw CSVs
├── rust/                   # Secondary Rust port (zero external crates)
│   ├── Cargo.toml
│   └── src/main.rs
├── results/
│   ├── environment.txt     # CPU model, OS, compiler version, optimization flags
│   └── merkle_results.csv  # Raw empirical benchmark measurements (1515 rows)
└── figures/                # Publication-quality 300 DPI evaluation plots
    ├── fig1_build_time_vs_K.png
    ├── fig2_proof_size_vs_N.png
    ├── fig3_adversarial_depth.png
    ├── fig4_multiproof_dedup.png
    ├── fig5_dmax_sweep.png
    └── fig6_pipeline_tarjan_knapsack.png
```

---

## Quickstart & Reproduction

### Prerequisites
- **C++ Compiler**: GCC 11+ or Clang 13+ or MSVC 2022 supporting `-std=c++20`.
- **Python**: Python 3.10+ with `pandas`, `numpy`, and `matplotlib`.
- **Rust (Optional)**: `cargo` 1.70+. If not installed, see [Installing Rust](#installing-rust) below.

### One-Command Run (Windows PowerShell)

```powershell
# 1. Compile C++ code
cd M:\DAA\merkle\cpp
.\build.bat

# 2. Run self-tests (validates SHA-256, Garsia-Wachs == Knuth, proofs, Tarjan, etc.)
.\build\merkle_tool.exe --selftest

# 3. Run all experiments E1 through E9 (writes to results/merkle_results.csv)
.\build\merkle_tool.exe --bench M:\DAA\merkle\results\merkle_results.csv

# 4. Generate all figures and summary tables
cd ..
& "C:\Users\Sapra\AppData\Local\Programs\Python\Python313\python.exe" python/plots.py
& "C:\Users\Sapra\AppData\Local\Programs\Python\Python313\python.exe" python/make_tables.py
```

*Expected Execution Time*: **~40 seconds total**.

### One-Command Run (Linux / macOS)

```bash
cd merkle/cpp
chmod +x build.sh
./build.sh
./build/merkle_tool --selftest
./build/merkle_tool --bench ../results/merkle_results.csv
cd ..
python3 python/plots.py
python3 python/make_tables.py
```

---

## Installing Rust

The primary experiments and benchmarks run entirely in standalone C++20 without dependencies. A clean Rust port is provided in `rust/`. If `cargo` is missing on your machine, install it via `rustup`:

- **Windows**:
  ```powershell
  winget install --id Rustlang.Rustup
  # or download from https://rustup.rs
  ```
- **Linux / macOS**:
  ```bash
  curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs | sh
  source $HOME/.cargo/env
  ```

Once installed, build and test the Rust port:
```bash
cd merkle/rust
cargo run --release -- --test
cargo run --release -- --bench ../results/merkle_rust_results.csv
```

---

## Key Results Summary

- **Garsia-Wachs Speedup**: Produces tree costs that match Knuth interval DP **exactly down to the integer unit**, while running **$62.3\times$ faster** at $K=1024$ (0.69 ms vs 43.36 ms).
- **Proof Size Reduction**: Hot account proof sizes decrease from 20 hashes ($675$ bytes) down to **8.88 hashes** ($318$ bytes), achieving a **52.9% bandwidth reduction**.
- **Multiproof Compression**: Shared-hash deduplication removes **79.9%** of witness bytes for a 500-transaction block.
- **Pre-Block Pipeline**: Tarjan's SCC identifies wash-trade circular transfer rings in $<125$ µs on $1000$-vertex transfer graphs, while 2D Knapsack DP extracts $8.4\% - 14.6\%$ higher fee revenue than greedy fee-per-gas under gas constraints.
