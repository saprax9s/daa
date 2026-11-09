"""
Plot generation script for Merkle Tree benchmark results.
Reads results/merkle_results.csv and outputs figures to figures/*.png.
"""

import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    csv_path = os.path.join(project_root, "results", "merkle_results.csv")
    fig_dir = os.path.join(project_root, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} does not exist.")
        return

    df = pd.read_csv(csv_path)
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

    # =========================================================================
    # Figure 1: Build Time vs K (E1)
    # =========================================================================
    e1_time = df[(df['experiment'] == 'E1') & (df['metric'] == 'build_time_us')]
    if not e1_time.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        for algo, group in e1_time.groupby('algo'):
            meds = group.groupby('K')['value'].median()
            ax.plot(meds.index, meds.values, marker='o', linewidth=2, label=algo)
        ax.set_xscale('log', base=2)
        ax.set_yscale('log')
        ax.set_xlabel('Hot Set Size K (Leaves = 2K+1)', fontsize=12)
        ax.set_ylabel('Median Build Time (µs)', fontsize=12)
        ax.set_title('E1: Merkle Tree Construction Time vs Hot Set Size K', fontsize=13, fontweight='bold')
        ax.legend(title='Algorithm', frameon=True)
        ax.grid(True, which="both", ls="--", alpha=0.5)
        plt.tight_layout()
        fig_path = os.path.join(fig_dir, "fig1_build_time_vs_K.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"Saved: {fig_path}")

    # =========================================================================
    # Figure 2: Proof Size vs N (E2)
    # =========================================================================
    e2 = df[df['experiment'] == 'E2']
    if not e2.empty:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        e2_hashes = e2[e2['metric'] == 'proof_hashes']
        e2_bytes = e2[e2['metric'] == 'proof_bytes']

        algos = e2_hashes['algo'].unique()
        n_vals = e2_hashes['N'].unique()
        x = np.arange(len(n_vals))
        width = 0.2

        for i, algo in enumerate(algos):
            vals_h = [e2_hashes[(e2_hashes['algo'] == algo) & (e2_hashes['N'] == n)]['value'].mean() for n in n_vals]
            vals_b = [e2_bytes[(e2_bytes['algo'] == algo) & (e2_bytes['N'] == n)]['value'].mean() for n in n_vals]
            ax1.bar(x + i * width, vals_h, width, label=algo)
            ax2.bar(x + i * width, vals_b, width, label=algo)

        ax1.set_xticks(x + width * (len(algos) - 1) / 2)
        ax1.set_xticklabels([f"N={int(n):,}" for n in n_vals])
        ax1.set_ylabel('Proof Size (Hashes)', fontsize=12)
        ax1.set_title('Proof Size in Hashes (K=128)', fontsize=12, fontweight='bold')
        ax1.legend(frameon=True)
        ax1.grid(True, axis='y', ls="--", alpha=0.5)

        ax2.set_xticks(x + width * (len(algos) - 1) / 2)
        ax2.set_xticklabels([f"N={int(n):,}" for n in n_vals])
        ax2.set_ylabel('Proof Size (Bytes)', fontsize=12)
        ax2.set_title('Proof Size in Bytes (K=128)', fontsize=12, fontweight='bold')
        ax2.legend(frameon=True)
        ax2.grid(True, axis='y', ls="--", alpha=0.5)

        plt.suptitle('E2: Membership Proof Size: Interleaved vs Full Balanced Baseline (B1)', fontsize=13, fontweight='bold')
        plt.tight_layout()
        fig_path = os.path.join(fig_dir, "fig2_proof_size_vs_N.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"Saved: {fig_path}")

    # =========================================================================
    # Figure 3: Adversarial Depth Analysis (E3)
    # =========================================================================
    e3_budget = df[(df['experiment'] == 'E3') & (df['algo'] == 'attacker_budget') & (df['metric'] == 'max_depth')]
    if not e3_budget.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.plot(e3_budget['seed'], e3_budget['value'], marker='s', color='#d62728', linewidth=2, label='Max Tree Depth')
        ax.set_xlabel('Attacker Gas/Touch Budget', fontsize=12)
        ax.set_ylabel('Maximum Leaf Depth', fontsize=12)
        ax.set_title('E3: Maximum Tree Depth Reached vs Attacker Budget (N=64)', fontsize=13, fontweight='bold')
        ax.grid(True, ls="--", alpha=0.5)
        ax.legend(frameon=True)
        plt.tight_layout()
        fig_path = os.path.join(fig_dir, "fig3_adversarial_depth.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"Saved: {fig_path}")

    # =========================================================================
    # Figure 4: Multiproof Batch Size vs Dedup Ratio (E4)
    # =========================================================================
    e4 = df[df['experiment'] == 'E4']
    if not e4.empty:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        mp_h = e4[(e4['algo'] == 'multiproof') & (e4['metric'] == 'multiproof_hashes')]
        sp_h = e4[(e4['algo'] == 'separate_proofs') & (e4['metric'] == 'separate_hashes')]
        dedup = e4[(e4['algo'] == 'multiproof') & (e4['metric'] == 'deduplication_ratio')]

        ax1.plot(mp_h['m'], mp_h['value'], marker='o', color='#1f77b4', linewidth=2, label='Multiproof (Deduped)')
        ax1.plot(sp_h['m'], sp_h['value'], marker='x', color='#ff7f0e', linewidth=2, linestyle='--', label='Separate Proofs')
        ax1.set_xlabel('Batch Size m (Transactions)', fontsize=12)
        ax1.set_ylabel('Total Hashes in Witness', fontsize=12)
        ax1.set_title('Witness Hashes vs Batch Size m', fontsize=12, fontweight='bold')
        ax1.legend(frameon=True)
        ax1.grid(True, ls="--", alpha=0.5)

        ax2.plot(dedup['m'], dedup['value'] * 100.0, marker='d', color='#2ca02c', linewidth=2)
        ax2.set_xlabel('Batch Size m (Transactions)', fontsize=12)
        ax2.set_ylabel('Deduplication Savings (%)', fontsize=12)
        ax2.set_title('Hash Deduplication Efficiency', fontsize=12, fontweight='bold')
        ax2.grid(True, ls="--", alpha=0.5)

        plt.suptitle('E4: Batch Multiproof Verification Scaling', fontsize=13, fontweight='bold')
        plt.tight_layout()
        fig_path = os.path.join(fig_dir, "fig4_multiproof_dedup.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"Saved: {fig_path}")

    # =========================================================================
    # Figure 5: Depth-Capped DP Sweep (E6)
    # =========================================================================
    e6 = df[df['experiment'] == 'E6']
    if not e6.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        e6_avg = e6[e6['metric'] == 'avg_depth']
        ax.plot(e6_avg['dmax'], e6_avg['value'], marker='o', color='#9467bd', linewidth=2, label='Traffic-Weighted Avg Depth')
        ax.axvline(x=8, color='gray', linestyle=':', label=r'Theoretical Minimum Cap $\lceil \log_2(N) \rceil$')
        ax.set_xlabel('Depth Cap Bound D_max', fontsize=12)
        ax.set_ylabel('Average Depth (Weighted)', fontsize=12)
        ax.set_title('E6: Average Depth Trade-Off vs Worst-Case Cap Bound D_max', fontsize=13, fontweight='bold')
        ax.legend(frameon=True)
        ax.grid(True, ls="--", alpha=0.5)
        plt.tight_layout()
        fig_path = os.path.join(fig_dir, "fig5_dmax_sweep.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"Saved: {fig_path}")

    # =========================================================================
    # Figure 6: Pre-Block Pipeline (E9)
    # =========================================================================
    e9 = df[df['experiment'] == 'E9']
    if not e9.empty:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        scc = e9[(e9['algo'] == 'tarjan_scc') & (e9['metric'] == 'scc_time_us')]
        if not scc.empty:
            med_scc = scc.groupby('K')['value'].median()
            ax1.plot(med_scc.index, med_scc.values, marker='o', color='#8c564b', linewidth=2)
            ax1.set_xlabel('Vertices V (Account Transfer Graph)', fontsize=12)
            ax1.set_ylabel('Tarjan SCC Runtime (µs)', fontsize=12)
            ax1.set_title('Tarjan SCC Wash-Trade Detection Runtime', fontsize=12, fontweight='bold')
            ax1.grid(True, ls="--", alpha=0.5)

        gain = e9[(e9['algo'] == 'dp_gain') & (e9['metric'] == 'gain_pct')]
        if not gain.empty:
            ax2.bar([str(int(g)) for g in gain['dmax']], gain['value'], color='#17becf', width=0.4)
            ax2.set_xlabel('Block Gas Limit G', fontsize=12)
            ax2.set_ylabel('Fee Gain Over Greedy (%)', fontsize=12)
            ax2.set_title('Knapsack DP Gain vs Greedy Fee-Per-Gas', fontsize=12, fontweight='bold')
            ax2.grid(True, axis='y', ls="--", alpha=0.5)

        plt.suptitle('E9: Pre-Block Validation and Transaction Selection Pipeline', fontsize=13, fontweight='bold')
        plt.tight_layout()
        fig_path = os.path.join(fig_dir, "fig6_pipeline_tarjan_knapsack.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"Saved: {fig_path}")

    print("All Merkle figures generated successfully.")

if __name__ == "__main__":
    main()

