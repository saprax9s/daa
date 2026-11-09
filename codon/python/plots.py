"""
Plot generation script for Codon Design benchmark results.
Reads results/codon_results.csv and outputs figures to figures/*.png.
"""

import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    csv_path = os.path.join(project_root, "results", "codon_results.csv")
    fig_dir = os.path.join(project_root, "figures")
    os.makedirs(fig_dir, exist_ok=True)

    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} does not exist.")
        return

    df = pd.read_csv(csv_path)
    plt.style.use('seaborn-v0_8-whitegrid' if 'seaborn-v0_8-whitegrid' in plt.style.available else 'default')

    # =========================================================================
    # Figure 1: Optimality Gap and Validity Rate (E1)
    # =========================================================================
    e1 = df[df['experiment'] == 'E1']
    if not e1.empty:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))
        feas = e1[e1['metric'] == 'feasible']
        gap = e1[e1['metric'] == 'optimality_gap_pct']

        # Validity rate bar chart
        methods = ['exact_dp', 'greedy_m3', 'sample_reject_m4']
        method_labels = ['Exact DP (M0)', 'Greedy (M3)', 'Sample & Reject (M4)']
        rates = [feas[feas['method'] == m]['value'].mean() * 100.0 for m in methods]
        colors = ['#2ca02c', '#1f77b4', '#ff7f0e']

        ax1.bar(method_labels, rates, color=colors, width=0.5)
        ax1.set_ylabel('Validity / Feasibility Rate (%)', fontsize=12)
        ax1.set_title('Feasibility Success Rate (20 Seeded Instances)', fontsize=12, fontweight='bold')
        ax1.set_ylim(0, 110)
        ax1.grid(True, axis='y', ls="--", alpha=0.5)

        # Optimality gap boxplot
        gap_data = [gap[gap['method'] == m]['value'].dropna().values for m in ['greedy_m3', 'sample_reject_m4']]
        ax2.boxplot(gap_data, tick_labels=['Greedy (M3)', 'Sample & Reject (M4)'], patch_artist=True)
        ax2.set_ylabel('Optimality Gap vs DP (%)', fontsize=12)
        ax2.set_title('Optimality Gap Distribution vs Exact Optimum', fontsize=12, fontweight='bold')
        ax2.grid(True, axis='y', ls="--", alpha=0.5)

        plt.suptitle('E1: Exact DP vs Heuristic Baselines (n=15)', fontsize=13, fontweight='bold')
        plt.tight_layout()
        fig_path = os.path.join(fig_dir, "fig1_optimality_gap.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"Saved: {fig_path}")

    # =========================================================================
    # Figure 2: Runtime vs n (E2 Log-Log Plot)
    # =========================================================================
    e2 = df[df['experiment'] == 'E2']
    if not e2.empty:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))
        e2_time = e2[e2['metric'] == 'runtime_us']
        e2_states = e2[e2['metric'] == 'dp_states']

        for mset, grp in e2_time.groupby('motif_set'):
            meds = grp.groupby('n')['value'].median()
            ax1.plot(meds.index, meds.values, marker='o', linewidth=2, label=mset)
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlabel('Protein Length n (Residues)', fontsize=12)
        ax1.set_ylabel('Median DP Runtime (µs)', fontsize=12)
        ax1.set_title('DP Runtime vs Sequence Length n (Log-Log)', fontsize=12, fontweight='bold')
        ax1.legend(title='Motif Set', frameon=True)
        ax1.grid(True, which="both", ls="--", alpha=0.5)

        for mset, grp in e2_states.groupby('motif_set'):
            st = grp.groupby('n')['value'].first()
            ax2.plot(st.index, st.values, marker='s', linewidth=2, label=mset)
        ax2.set_xscale('log')
        ax2.set_yscale('log')
        ax2.set_xlabel('Protein Length n (Residues)', fontsize=12)
        ax2.set_ylabel('Total DP States Explored', fontsize=12)
        ax2.set_title('DP States vs Sequence Length n (Log-Log)', fontsize=12, fontweight='bold')
        ax2.legend(title='Motif Set', frameon=True)
        ax2.grid(True, which="both", ls="--", alpha=0.5)

        plt.suptitle('E2: Empirical Complexity vs Sequence Length (Theory: O(n² |Q|))', fontsize=13, fontweight='bold')
        plt.tight_layout()
        fig_path = os.path.join(fig_dir, "fig2_runtime_vs_n.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"Saved: {fig_path}")

    # =========================================================================
    # Figure 3: Runtime and States vs |Q| (E3)
    # =========================================================================
    e3 = df[df['experiment'] == 'E3']
    if not e3.empty:
        fig, ax1 = plt.subplots(figsize=(8, 5))
        e3_time = e3[e3['metric'] == 'runtime_us']
        e3_states = e3[e3['metric'] == 'automaton_states']
        
        q_vals = []
        runtimes = []
        dp_states = []
        for mset, grp in e3.groupby('motif_set'):
            q_row = grp[grp['metric'] == 'automaton_states']
            t_row = grp[grp['metric'] == 'runtime_us']
            s_row = grp[grp['metric'] == 'dp_states']
            if not q_row.empty and not t_row.empty:
                q_vals.append(q_row['value'].values[0])
                runtimes.append(t_row['value'].median())
                dp_states.append(s_row['value'].values[0] if not s_row.empty else 0)

        # Sort by |Q|
        order = np.argsort(q_vals)
        q_vals = np.array(q_vals)[order]
        runtimes = np.array(runtimes)[order]

        ax1.plot(q_vals, runtimes, marker='o', color='#1f77b4', linewidth=2, label='Runtime (µs)')
        ax1.set_xlabel('Automaton State Count |Q|', fontsize=12)
        ax1.set_ylabel('DP Runtime (µs)', fontsize=12, color='#1f77b4')
        ax1.set_title('E3: DP Runtime vs Automaton Size |Q| at Fixed n=50', fontsize=13, fontweight='bold')
        ax1.grid(True, ls="--", alpha=0.5)
        plt.tight_layout()
        fig_path = os.path.join(fig_dir, "fig3_runtime_vs_Q.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"Saved: {fig_path}")

    # =========================================================================
    # Figure 4: Pruning Ablation (E4)
    # =========================================================================
    e4 = df[df['experiment'] == 'E4']
    if not e4.empty:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))
        modes = ['D0_no_pruning', 'D1_gc_pruning', 'D2_dominance']
        labels = ['D0 (None)', 'D1 (GC Bound)', 'D2 (Dominance)']
        
        times = [e4[(e4['method'] == m) & (e4['metric'] == 'runtime_us')]['value'].median() for m in modes]
        states = [e4[(e4['method'] == m) & (e4['metric'] == 'dp_states')]['value'].iloc[0] for m in modes]

        colors = ['#d62728', '#2ca02c', '#1f77b4']
        ax1.bar(labels, states, color=colors, width=0.5)
        ax1.set_ylabel('Total DP States Explored', fontsize=12)
        ax1.set_title('State Space Reduction (n=30)', fontsize=12, fontweight='bold')
        ax1.grid(True, axis='y', ls="--", alpha=0.5)

        ax2.bar(labels, times, color=colors, width=0.5)
        ax2.set_ylabel('Median Runtime (µs)', fontsize=12)
        ax2.set_title('Execution Time Impact (n=30)', fontsize=12, fontweight='bold')
        ax2.grid(True, axis='y', ls="--", alpha=0.5)

        plt.suptitle('E4: Pruning Ablation on D0 vs D1 vs D2', fontsize=13, fontweight='bold')
        plt.tight_layout()
        fig_path = os.path.join(fig_dir, "fig4_pruning_ablation.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"Saved: {fig_path}")

    # =========================================================================
    # Figure 5: Window GC Scaling (E5)
    # =========================================================================
    e5 = df[df['experiment'] == 'E5']
    if not e5.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        k_vals = []
        states_vals = []
        times_vals = []
        for method, grp in e5.groupby('method'):
            # method is e.g. window_gc_k3
            k = int(method.replace('window_gc_k', ''))
            st = grp[grp['metric'] == 'dp_states']['value'].values[0]
            tm = grp[grp['metric'] == 'runtime_us']['value'].median()
            k_vals.append(k)
            states_vals.append(st)
            times_vals.append(tm)
            
        order = np.argsort(k_vals)
        k_vals = np.array(k_vals)[order]
        states_vals = np.array(states_vals)[order]
        times_vals = np.array(times_vals)[order]

        ax.plot(k_vals, states_vals, marker='o', color='#e377c2', linewidth=2, label='DP States')
        ax.set_xlabel('Window Size k (Codons)', fontsize=12)
        ax.set_ylabel('Total DP States Explored', fontsize=12)
        ax.set_title('E5: Window GC State Growth vs Window Size k (n=20)', fontsize=13, fontweight='bold')
        ax.grid(True, ls="--", alpha=0.5)
        ax.legend(frameon=True)
        plt.tight_layout()
        fig_path = os.path.join(fig_dir, "fig5_window_gc.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"Saved: {fig_path}")

    # =========================================================================
    # Figure 6: Realistic Proteins GFP & Insulin (E6)
    # =========================================================================
    e6 = df[df['experiment'] == 'E6']
    if not e6.empty:
        fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11, 5))
        e6_time = e6[e6['metric'] == 'runtime_us']
        e6_states = e6[e6['metric'] == 'dp_states']

        labels = []
        t_vals = []
        s_vals = []
        for (meth, gc), grp in e6_time.groupby(['method', 'gc_mode']):
            labels.append(f"{meth}\n({gc})")
            t_vals.append(grp['value'].median())
            st_row = e6_states[(e6_states['method'] == meth) & (e6_states['gc_mode'] == gc)]
            s_vals.append(st_row['value'].values[0] if not st_row.empty else 0)

        x = np.arange(len(labels))
        ax1.bar(x, t_vals, color='#17becf', width=0.5)
        ax1.set_xticks(x)
        ax1.set_xticklabels(labels, rotation=15)
        ax1.set_ylabel('Median Runtime (µs)', fontsize=12)
        ax1.set_title('Design Runtime under 15 Restriction Sites', fontsize=12, fontweight='bold')
        ax1.grid(True, axis='y', ls="--", alpha=0.5)

        ax2.bar(x, s_vals, color='#bcbd22', width=0.5)
        ax2.set_xticks(x)
        ax2.set_xticklabels(labels, rotation=15)
        ax2.set_ylabel('DP States Explored', fontsize=12)
        ax2.set_title('DP States under 15 Restriction Sites', fontsize=12, fontweight='bold')
        ax2.grid(True, axis='y', ls="--", alpha=0.5)

        plt.suptitle('E6: Realistic Restriction-Site Avoidance on Insulin (51 aa) and GFP (238 aa)', fontsize=13, fontweight='bold')
        plt.tight_layout()
        fig_path = os.path.join(fig_dir, "fig6_realistic_proteins.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"Saved: {fig_path}")

    # =========================================================================
    # Figure 7: Sample-and-Reject Wasted Candidates (E7)
    # =========================================================================
    e7 = df[df['experiment'] == 'E7']
    if not e7.empty:
        fig, ax = plt.subplots(figsize=(8, 5))
        e7_wasted = e7[e7['metric'] == 'wasted_pct']
        margins = []
        w_pct = []
        for gc, grp in e7_wasted.groupby('gc_mode'):
            m = int(gc.replace('margin_', ''))
            margins.append(m)
            w_pct.append(grp['value'].values[0])
        order = np.argsort(margins)
        margins = np.array(margins)[order]
        w_pct = np.array(w_pct)[order]

        ax.plot(margins * 2, w_pct, marker='s', color='#d62728', linewidth=2)
        ax.set_xlabel('GC Window Allowed Range Width (Nucleotides)', fontsize=12)
        ax.set_ylabel('Wasted Candidates (%)', fontsize=12)
        ax.set_title('E7: Sample-and-Reject Rejection Rate vs Constraint Tightness', fontsize=13, fontweight='bold')
        ax.grid(True, ls="--", alpha=0.5)
        plt.tight_layout()
        fig_path = os.path.join(fig_dir, "fig7_sample_wasted.png")
        plt.savefig(fig_path, dpi=300)
        plt.close()
        print(f"Saved: {fig_path}")

    print("All Codon figures generated successfully.")

if __name__ == "__main__":
    main()
