"""
Table generation script for Codon Design benchmark results.
Reads results/codon_results.csv and formats summary markdown tables for docs/REPORT.md.
"""

import os
import pandas as pd
import numpy as np

def to_md(df):
    headers = list(df.columns)
    lines = ["| " + " | ".join(str(h) for h in headers) + " |"]
    lines.append("| " + " | ".join(["---"] * len(headers)) + " |")
    for _, row in df.iterrows():
        lines.append("| " + " | ".join(str(val) for val in row) + " |")
    return "\n".join(lines)

def format_e1_table(df):
    e1 = df[df['experiment'] == 'E1']
    if e1.empty:
        return ""
    feas = e1[e1['metric'] == 'feasible']
    gap = e1[e1['metric'] == 'optimality_gap_pct']
    records = []
    for method in ['exact_dp', 'greedy_m3', 'sample_reject_m4']:
        f_sub = feas[feas['method'] == method]
        g_sub = gap[gap['method'] == method]
        rate = f_sub['value'].mean() * 100.0 if not f_sub.empty else 0.0
        med_gap = g_sub['value'].median() if not g_sub.empty else 0.0
        mean_gap = g_sub['value'].mean() if not g_sub.empty else 0.0
        records.append({
            'Method': method,
            'Instances': len(f_sub),
            'Feasibility Rate (%)': f"{rate:.1f}%",
            'Median Gap vs DP (%)': f"{med_gap:.2f}%",
            'Mean Gap vs DP (%)': f"{mean_gap:.2f}%"
        })
    return to_md(pd.DataFrame(records))

def format_e2_table(df):
    e2 = df[df['experiment'] == 'E2']
    if e2.empty:
        return ""
    e2_time = e2[e2['metric'] == 'runtime_us']
    records = []
    for (mset, n), grp in e2_time.groupby(['motif_set', 'n']):
        vals = grp['value']
        med = np.median(vals)
        mn = np.min(vals)
        iqr = np.percentile(vals, 75) - np.percentile(vals, 25)
        st_row = df[(df['experiment'] == 'E2') & (df['motif_set'] == mset) & (df['n'] == n) & (df['metric'] == 'dp_states')]
        pk_row = df[(df['experiment'] == 'E2') & (df['motif_set'] == mset) & (df['n'] == n) & (df['metric'] == 'peak_states')]
        states = st_row['value'].values[0] if not st_row.empty else 0
        peak = pk_row['value'].values[0] if not pk_row.empty else 0
        records.append({
            'Motif Set': mset,
            'n (Residues)': int(n),
            'Median Runtime (µs)': f"{med:.2f}",
            'Min (µs)': f"{mn:.2f}",
            'IQR (µs)': f"{iqr:.2f}",
            'Total DP States': f"{int(states):,}",
            'Peak Layer States': f"{int(peak):,}"
        })
    return to_md(pd.DataFrame(records))

def format_e4_table(df):
    e4 = df[df['experiment'] == 'E4']
    if e4.empty:
        return ""
    records = []
    for method, grp in e4.groupby('method'):
        t_row = grp[grp['metric'] == 'runtime_us']
        s_row = grp[grp['metric'] == 'dp_states']
        sc_row = grp[grp['metric'] == 'score']
        med_t = t_row['value'].median() if not t_row.empty else 0.0
        states = s_row['value'].values[0] if not s_row.empty else 0
        score = sc_row['value'].values[0] if not sc_row.empty else 0
        records.append({
            'Configuration': method,
            'Total States': f"{int(states):,}",
            'Median Runtime (µs)': f"{med_t:.2f}",
            'Optimal Score': int(score)
        })
    return to_md(pd.DataFrame(records))

def format_e6_table(df):
    e6 = df[df['experiment'] == 'E6']
    if e6.empty:
        return ""
    records = []
    for (method, gc), grp in e6.groupby(['method', 'gc_mode']):
        t_row = grp[grp['metric'] == 'runtime_us']
        f_row = grp[grp['metric'] == 'feasible']
        sc_row = grp[grp['metric'] == 'score']
        s_row = grp[grp['metric'] == 'dp_states']
        med_t = t_row['value'].median() if not t_row.empty else 0.0
        feas = bool(f_row['value'].values[0]) if not f_row.empty else False
        score = sc_row['value'].values[0] if not sc_row.empty else 0
        states = s_row['value'].values[0] if not s_row.empty else 0
        records.append({
            'Protein Target': method,
            'GC Constraint': gc,
            'Feasible': "Yes" if feas else "No",
            'Optimal Score': f"{int(score):,}",
            'States Explored': f"{int(states):,}",
            'Runtime (µs)': f"{med_t:.2f}"
        })
    return to_md(pd.DataFrame(records))

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    csv_path = os.path.join(project_root, "results", "codon_results.csv")

    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} not found.")
        return

    df = pd.read_csv(csv_path)
    print("### Table 1: Optimality Gap and Validity Rate (E1)\n")
    print(format_e1_table(df))
    print("\n### Table 2: Scaling vs Sequence Length n (E2)\n")
    print(format_e2_table(df))
    print("\n### Table 3: Pruning Ablation D0 / D1 / D2 (E4)\n")
    print(format_e4_table(df))
    print("\n### Table 4: Realistic Protein Design: GFP and Insulin (E6)\n")
    print(format_e6_table(df))

if __name__ == "__main__":
    main()

