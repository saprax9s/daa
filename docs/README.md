# DAA Case Studies — Repository Index

Comprehensive implementation, rigorous empirical benchmarking, and documentation of two case studies for the Design and Analysis of Algorithms (DAA) course (single author; no new algorithm claimed; rigorous empirical investigation).

---

## 1. Project 1: Exact Motif-Free Codon Design (`codon/`)

- **Topic**: Exact polynomial-time reverse-translation of protein sequences avoiding all forbidden restriction sites and motifs, satisfying global/sliding-window GC constraints, and maximizing an additive codon/pair score.
- **Specification**: [`codon/prompt.md`](file:///M:/DAA/codon/prompt.md)
- **Documentation**: [`codon/docs/REPORT.md`](file:///M:/DAA/codon/docs/REPORT.md)
- **Baselines**: [`codon/BASELINES.md`](file:///M:/DAA/codon/BASELINES.md)
- **Status**: **Complete**. Passed validations V1–V6. All experiments E1–E8 benchmarked (1,283 CSV rows), 7 publication figures generated in `codon/figures/`, and living report complete.

### Quick Run (`codon/`)
```powershell
cd M:\DAA\codon\cpp
.\build.bat
.\build\codon_tool.exe --selftest
.\build\codon_tool.exe --generate-instances M:\DAA\codon\instances
.\build\codon_tool.exe --bench M:\DAA\codon\results\codon_results.csv
cd ..
& "C:\Users\Sapra\AppData\Local\Programs\Python\Python313\python.exe" python/plots.py
& "C:\Users\Sapra\AppData\Local\Programs\Python\Python313\python.exe" python/make_tables.py
```

---

## 2. Project 2: Depth-Bounded Alphabetic Merkle Trees (`merkle/`)

- **Topic**: Depth-bounded alphabetic Merkle trees for stateless blockchain verification under Zipf-skewed account traffic ($s \ge 1.0$), comparing 6 tree structures with a hard depth ceiling $D_{\max}$ and pre-block pipeline validation (Tarjan SCC & 2D Knapsack DP).
- **Specification**: [`merkle/prompt.md`](file:///M:/DAA/merkle/prompt.md)
- **Documentation**: [`merkle/docs/REPORT.md`](file:///M:/DAA/merkle/docs/REPORT.md)
- **Baselines**: [`merkle/BASELINES.md`](file:///M:/DAA/merkle/BASELINES.md)
- **Status**: **Complete**. Passed all 17 selftests. All experiments E1–E9 benchmarked (1,515 CSV rows), 6 publication figures generated in `merkle/figures/`, and living report complete.

### Quick Run (`merkle/`)
```powershell
cd M:\DAA\merkle\cpp
.\build.bat
.\build\merkle_tool.exe --selftest
.\build\merkle_tool.exe --bench M:\DAA\merkle\results\merkle_results.csv
cd ..
& "C:\Users\Sapra\AppData\Local\Programs\Python\Python313\python.exe" python/plots.py
& "C:\Users\Sapra\AppData\Local\Programs\Python\Python313\python.exe" python/make_tables.py
```

---

## Standard Repository Layout (Both Projects)

```
project/
  README.md            exact run instructions (Windows PowerShell + Linux)
  BASELINES.md         every baseline: what, why, parameters, fairness notes
  docs/REPORT.md       living empirical report populated with genuine data
  cpp/                 src/, tests/, bench/, CMakeLists.txt, build.bat, build.sh
  python/              oracle.py, plots.py, make_tables.py
  rust/                standalone Rust port (zero external dependencies)
  instances/           generated *.json test/benchmark instances
  results/             raw CSVs & environment.txt
  figures/             publication-quality 300 DPI plots
```

---

## Environment & Toolchain
- **CPU**: AMD Ryzen 7 7840HS w/ Radeon 780M Graphics (8 cores, 16 threads)
- **OS**: Windows 11 Home (64-bit)
- **Compiler**: MinGW-w64 GCC 15.2.0 (`-std=c++20 -O2 -Wall -Wextra`)
- **Python**: Python 3.13 (`pandas`, `numpy`, `matplotlib`)
