#!/usr/bin/env python3
"""
Analyze SASA differences between:
1. Modified Cysteines (Glut) that are CLOSE (<=10A) to a crosstalk PTM.
2. Modified Cysteines (Glut) that are FAR (>10A) from any crosstalk PTM.
3. Control Cysteines (Non-modified).

Inputs:
  --distances reports/distances.csv
  --sasa reports/sasa_cys.csv
Output:
  --out-csv reports/ptm_sasa_summary.csv
  --out-md reports/sasa_summary.md
"""
import argparse
import pandas as pd
import numpy as np
from pathlib import Path

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--distances", required=True)
    p.add_argument("--sasa", required=True)
    p.add_argument("--out-csv", default="reports/ptm_sasa_summary.csv")
    p.add_argument("--out-md", default="reports/sasa_summary.md")
    args = p.parse_args()

    # Load Data
    if not Path(args.distances).exists():
        print("Distances file not found.")
        return
    if not Path(args.sasa).exists():
        print("SASA file not found.")
        return

    df_dist = pd.read_csv(args.distances)
    df_sasa = pd.read_csv(args.sasa)

    # Ensure types
    df_sasa["UniProt_ID"] = df_sasa["UniProt_ID"].astype(str).str.strip()
    df_sasa["resnum"] = df_sasa["resnum"].astype(int)
    
    df_dist["UniProt_ID"] = df_dist["UniProt_ID"].astype(str).str.strip()
    df_dist["cysteine_pos"] = df_dist["cysteine_pos"].astype(int)

    # Aggregating distances per cysteine
    # We want the MINIMUM distance to ANY PTM for each cysteine
    # This determines if it's "Close" or "Far" from the PTM environment
    
    # Group by Cysteine (ID + Pos)
    # Calculate min distance for this cysteine across all pairs found
    cys_groups = df_dist.groupby(["UniProt_ID", "cysteine_pos", "label"])
    
    cys_summary = cys_groups.agg(
        min_dist_to_ptm=("min_distance_A", "min"),
        pair_count=("ptm_type", "count")
    ).reset_index()

    # Merge with SASA
    # Left join: We only care about cysteines that were analyzed for distances
    merged = pd.merge(
        cys_summary,
        df_sasa[["UniProt_ID", "resnum", "SASA"]],
        left_on=["UniProt_ID", "cysteine_pos"],
        right_on=["UniProt_ID", "resnum"],
        how="left"
    )

    # Categorize
    def categorize(row):
        if row["label"] == "control":
            return "Control (Non-modified)"
        # It is Glut
        if row["min_dist_to_ptm"] <= 10.0:
            return "Glut (Close to PTM)"
        else:
            return "Glut (Far from PTM)"

    merged["Group"] = merged.apply(categorize, axis=1)

    # Statistics
    stats = merged.groupby("Group")["SASA"].agg(["count", "mean", "median", "std"]).reset_index()
    stats = stats.sort_values("mean", ascending=False)

    # Write CSV
    merged.to_csv(args.out_csv, index=False)

    # Write Markdown Report
    md_lines = []
    md_lines.append("### SASA Analysis Summary")
    md_lines.append("")
    md_lines.append("| Group | Count | Mean SASA (Å²) | Median SASA (Å²) | Std Dev |")
    md_lines.append("|---|---|---|---|---|")
    
    for _, row in stats.iterrows():
        md_lines.append(f"| {row['Group']} | {int(row['count'])} | {row['mean']:.2f} | {row['median']:.2f} | {row['std']:.2f} |")
    
    md_lines.append("")
    md_lines.append("**Definitions:**")
    md_lines.append("- **Glut (Close to PTM)**: Modified cysteines within 10Å of any crosstalk PTM.")
    md_lines.append("- **Glut (Far from PTM)**: Modified cysteines farther than 10Å from all crosstalk PTMs.")
    md_lines.append("- **Control**: Non-modified cysteines on the same proteins.")
    
    with open(args.out_md, "w") as f:
        f.write("\n".join(md_lines))

    print(f"SASA analysis written to {args.out_csv} and {args.out_md}")

if __name__ == "__main__":
    main()
