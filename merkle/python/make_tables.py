"""
Table generation script for Merkle Tree benchmark results.
Reads results/merkle_results.csv and formats summary markdown tables for docs/REPORT.md.
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
    e1 = df[(df['experiment'] == 'E1') & (df['metric'] == 'build_time_us')]
    if e1.empty:
        return ""
    
    records = []
    for (algo, k), grp in e1.groupby(['algo', 'K']):
        vals = grp['value']
        med = np.median(vals)
        mn = np.min(vals)
        iqr = np.percentile(vals, 75) - np.percentile(vals, 25)
        cost_row = df[(df['experiment'] == 'E1') & (df['algo'] == algo) & (df['K'] == k) & (df['metric'] == 'cost')]
        depth_row = df[(df['experiment'] == 'E1') & (df['algo'] == algo) & (df['K'] == k) & (df['metric'] == 'max_depth')]
        cost = cost_row['value'].values[0] if not cost_row.empty else 0
        depth = depth_row['value'].values[0] if not depth_row.empty else 0
        records.append({
            'Algorithm': algo,
            'K': int(k),
            'Leaves': int(2 * k + 1),
            'Median (µs)': f"{med:.2f}",
            'Min (µs)': f"{mn:.2f}",
            'IQR (µs)': f"{iqr:.2f}",
            'Cost': f"{int(cost):,}",
            'Max Depth': int(depth)
        })
    res_df = pd.DataFrame(records)
    return to_md(res_df)

def format_e2_table(df):
    e2 = df[df['experiment'] == 'E2']
    if e2.empty:
        return ""
    records = []
    for (algo, n), grp in e2.groupby(['algo', 'N']):
        h_row = grp[grp['metric'] == 'proof_hashes']
        b_row = grp[grp['metric'] == 'proof_bytes']
        h_val = h_row['value'].values[0] if not h_row.empty else 0
        b_val = b_row['value'].values[0] if not b_row.empty else 0
        records.append({
            'Method': algo,
            'N (Accounts)': f"{int(n):,}",
            'Avg Hashes': f"{h_val:.2f}",
            'Avg Bytes': f"{b_val:.2f}"
        })
    return to_md(pd.DataFrame(records))

def format_e4_table(df):
    e4 = df[df['experiment'] == 'E4']
    if e4.empty:
        return ""
    records = []
    for m, grp in e4.groupby('m'):
        mp_h = grp[(grp['algo'] == 'multiproof') & (grp['metric'] == 'multiproof_hashes')]['value'].values[0]
        mp_b = grp[(grp['algo'] == 'multiproof') & (grp['metric'] == 'multiproof_bytes')]['value'].values[0]
        sp_h = grp[(grp['algo'] == 'separate_proofs') & (grp['metric'] == 'separate_hashes')]['value'].values[0]
        sp_b = grp[(grp['algo'] == 'separate_proofs') & (grp['metric'] == 'separate_bytes')]['value'].values[0]
        dedup = grp[(grp['algo'] == 'multiproof') & (grp['metric'] == 'deduplication_ratio')]['value'].values[0]
        records.append({
            'Batch m': int(m),
            'Multiproof Hashes': int(mp_h),
            'Separate Hashes': int(sp_h),
            'Multiproof Bytes': int(mp_b),
            'Separate Bytes': int(sp_b),
            'Savings (%)': f"{dedup * 100:.1f}%"
        })
    return to_md(pd.DataFrame(records))

def main():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(script_dir)
    csv_path = os.path.join(project_root, "results", "merkle_results.csv")

    if not os.path.exists(csv_path):
        print(f"Error: {csv_path} not found.")
        return

    df = pd.read_csv(csv_path)
    print("### Table 1: Construction Performance (E1)\n")
    print(format_e1_table(df))
    print("\n### Table 2: Membership Proof Sizes vs N (E2)\n")
    print(format_e2_table(df))
    print("\n### Table 3: Multiproof Batch Scaling (E4)\n")
    print(format_e4_table(df))

if __name__ == "__main__":
    main()
