#!/usr/bin/env python3
"""
Enrichment of other PTMs near glutathionylated cysteines vs unmodified cysteines.

Inputs
- reports/distances.csv (from scripts/compute_distances.py)
- Dataset/cleaned_glutathionylation_sites.csv (UniProt_ID, Glutathione_Site)
- alphafold_structures/AF-<UniProt>-F1-*.pdb (for control cysteine discovery)

Outputs
- reports/enrichment_summary.md: odds ratios and two-sided Fisher p-values
- reports/cysteine_summary.csv: per-cysteine proximity features
"""
from __future__ import annotations

import argparse
from pathlib import Path
import math
import pandas as pd
import numpy as np
import sys

# Local import path setup
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from scripts.pdb_utils import load_residue_letters, find_alphafold_pdb


def haldane_anscombe_or(a, b, c, d):
    a += 0.5; b += 0.5; c += 0.5; d += 0.5
    return (a / b) / (c / d)


def logC(n, k):
    from math import lgamma
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


def fisher_two_sided(a, b, c, d):
    # Based on hypergeometric; fixed margins
    n1 = a + b
    n2 = c + d
    m1 = a + c
    N = n1 + n2
    # probability of a given x in row1 successes
    def logp(x):
        return logC(m1, x) + logC(N - m1, n1 - x) - logC(N, n1)
    obs = a
    lp_obs = logp(obs)
    # enumerate feasible x
    lo = max(0, n1 - (N - m1))
    hi = min(n1, m1)
    s = 0.0
    for x in range(lo, hi + 1):
        lp = logp(x)
        if lp <= lp_obs + 1e-12:
            s += math.exp(lp)
    # numerical guard: cap at 1
    return min(1.0, s)


def build_control_cysteines(ids):
    controls = {}
    pdb_dir = Path('alphafold_structures')
    for uid in ids:
        pdb = find_alphafold_pdb(uid, pdb_dir)
        if not pdb or not pdb.exists():
            continue
        letters = load_residue_letters(pdb)
        c_positions = [pos for pos, aa in letters.items() if aa == 'C']
        controls[uid] = set(c_positions)
    return controls


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--dist', default='reports/distances.csv')
    ap.add_argument('--glut', default='Dataset/cleaned_glutathionylation_sites.csv')
    ap.add_argument('--out-md', default='reports/enrichment_summary.md')
    ap.add_argument('--out-csv', default='reports/cysteine_summary.csv')
    args = ap.parse_args()

    dists = pd.read_csv(args.dist)
    glut = pd.read_csv(args.glut)
    glut = glut[["UniProt_ID", "Glutathione_Site"]].copy()
    glut["Glutathione_Site"] = glut["Glutathione_Site"].astype(int)

    # Proteins in analysis
    proteins = sorted(set(dists["UniProt_ID"]))

    # Per-cysteine summary: aggregate min distances and presence flags by PTM type
    agg = dists.groupby(["UniProt_ID", "cysteine_pos", "ptm_type"]).agg(
        min_d=("min_distance_A", "min"),
        any5=("d_le_5A", "max"),
        any8=("d_le_8A", "max"),
        any10=("d_le_10A", "max"),
    ).reset_index()

    # Any PTM presence by cysteine
    # Include label if present
    if 'label' in dists.columns:
        labels = dists[['UniProt_ID','cysteine_pos','label']].drop_duplicates()
        agg = agg.merge(labels, on=['UniProt_ID','cysteine_pos'], how='left')
    any_by_cys = agg.groupby(["UniProt_ID", "cysteine_pos"]).agg(
        any5=("any5", "max"), any8=("any8", "max"), any10=("any10", "max")
    ).reset_index()
    if 'label' in dists.columns:
        any_by_cys = any_by_cys.merge(labels, on=['UniProt_ID','cysteine_pos'], how='left')

    # Build set of glut cysteines (only those with PDB present in dists)
    glut_set = set((r.UniProt_ID, int(r.Glutathione_Site)) for r in glut.itertuples(index=False))
    glut_set = {(u, p) for (u, p) in glut_set if u in proteins}
    if 'label' in dists.columns:
        # Derive control set directly from distances input
        all_pairs = set((r.UniProt_ID, int(r.cysteine_pos)) for r in dists[['UniProt_ID','cysteine_pos']].drop_duplicates().itertuples(index=False))
        control_set = {(u,p) for (u,p) in all_pairs if (u,p) not in glut_set}
    else:
        controls_map = build_control_cysteines(proteins)
        control_set = set()
        for u, poss in controls_map.items():
            for p in poss:
                if (u, p) not in glut_set:
                    control_set.add((u, p))

    # Merge presence flags for 'any PTM'
    any_lookup = {(r.UniProt_ID, int(r.cysteine_pos)): (int(r.any5), int(r.any8), int(r.any10))
                  for r in any_by_cys.itertuples(index=False)}

    def get_flags(pair):
        return any_lookup.get(pair, (0, 0, 0))

    # Build per-cysteine summary CSV rows (glut + controls)
    rows = []
    for pair in sorted(glut_set):
        a5, a8, a10 = get_flags(pair)
        rows.append({"UniProt_ID": pair[0], "cysteine_pos": pair[1], "label": "glut", "any5": a5, "any8": a8, "any10": a10})
    for pair in sorted(control_set):
        a5, a8, a10 = get_flags(pair)
        rows.append({"UniProt_ID": pair[0], "cysteine_pos": pair[1], "label": "control", "any5": a5, "any8": a8, "any10": a10})
    per_cys = pd.DataFrame(rows)
    Path(args.out_csv).parent.mkdir(parents=True, exist_ok=True)
    per_cys.to_csv(args.out_csv, index=False)

    # Enrichment calculations (any PTM and by PTM type)
    def compute_counts(per_cys_df, flag_col):
        g = per_cys_df[per_cys_df.label == 'glut']
        c = per_cys_df[per_cys_df.label == 'control']
        a = int(g[flag_col].sum())
        b = int((g.shape[0]) - a)
        c1 = int(c[flag_col].sum())
        d = int((c.shape[0]) - c1)
        return a, b, c1, d

    # Overall any-PTM
    results = []
    for flag in ["any5", "any8", "any10"]:
        a, b, c1, d = compute_counts(per_cys, flag)
        OR = haldane_anscombe_or(a, b, c1, d)
        p = fisher_two_sided(a, b, c1, d)
        results.append(("any", flag, a, b, c1, d, OR, p))

    # By PTM type presence
    for ptm in sorted(dists["ptm_type"].unique()):
        pres = agg[agg.ptm_type == ptm].groupby(["UniProt_ID", "cysteine_pos"]).agg(
            any5=("any5", "max"), any8=("any8", "max"), any10=("any10", "max")
        ).reset_index()
        pres_lookup = {(r.UniProt_ID, int(r.cysteine_pos)): (int(r.any5), int(r.any8), int(r.any10))
                       for r in pres.itertuples(index=False)}
        rows2 = []
        for pair in sorted(glut_set):
            a5, a8, a10 = pres_lookup.get(pair, (0, 0, 0))
            rows2.append({"label": "glut", "any5": a5, "any8": a8, "any10": a10})
        for pair in sorted(control_set):
            a5, a8, a10 = pres_lookup.get(pair, (0, 0, 0))
            rows2.append({"label": "control", "any5": a5, "any8": a8, "any10": a10})
        df2 = pd.DataFrame(rows2)
        for flag in ["any5", "any8", "any10"]:
            a, b, c1, d = compute_counts(df2, flag)
            OR = haldane_anscombe_or(a, b, c1, d)
            p = fisher_two_sided(a, b, c1, d)
            results.append((ptm, flag, a, b, c1, d, OR, p))

    # Write markdown summary
    lines = ["# PTM Proximity Enrichment", "", "Compares glutathionylated cysteines vs unmodified cysteines on the same proteins.", ""]
    lines.append("## Summary (odds ratio; Fisher two-sided p)")
    for ptm, flag, a, b, c1, d, OR, p in results:
        thr = {"any5": "<=5 Å", "any8": "<=8 Å", "any10": "<=10 Å"}[flag]
        lines.append(f"- {ptm} {thr}: OR={OR:.2f}, p={p:.2e} [glut {a}/{a+b}; ctrl {c1}/{c1+d}]")

    Path(args.out_md).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_md).write_text("\n".join(lines) + "\n")
    print(f"Wrote {args.out_md} and {args.out_csv}")


if __name__ == "__main__":
    raise SystemExit(main())
