#!/usr/bin/env python3
"""
Compute all-atom minimum distances between glutathionylated cysteines and
other PTM residues on the same protein using AlphaFold PDBs.

Inputs
- Dataset/cleaned_glutathionylation_sites.csv (UniProt_ID, Glutathione_Site)
- reports/artifacts/ptm_sites.filtered.csv (from scripts/ptm_loader.py)
- alphafold_structures/AF-<UniProt>-F1-*.pdb files

Output
- reports/distances.csv with columns:
  UniProt_ID, cysteine_pos, ptm_type, ptm_pos, min_distance_A
  plus flags for distance thresholds (<=5A, <=8A, <=10A).
"""
from __future__ import annotations

import argparse
from pathlib import Path
import math
import pandas as pd
import sys
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from scripts.pdb_utils import (
    load_residue_atoms,
    min_residue_distance,
    load_residue_letters,
    find_alphafold_pdb,
)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--glut", default="Dataset/cleaned_glutathionylation_sites.csv")
    ap.add_argument("--ptm", default="reports/artifacts/ptm_sites.filtered.csv")
    ap.add_argument("--pdb", default="alphafold_structures")
    ap.add_argument("--out", default="reports/distances.csv")
    ap.add_argument("--include-controls", action="store_true", help="Also compute distances for non-glut cysteines as controls")
    args = ap.parse_args()

    gout = pd.read_csv(args.glut)
    if not {"UniProt_ID", "Glutathione_Site"}.issubset(gout.columns):
        raise ValueError("Glutathionylation CSV must have UniProt_ID and Glutathione_Site")
    ptm = pd.read_csv(args.ptm)
    req_cols = {"UniProt_ID", "PTM", "Position"}
    if not req_cols.issubset(ptm.columns):
        raise ValueError(f"PTM CSV missing columns: {req_cols - set(ptm.columns)}")

    pdb_dir = Path(args.pdb)

    # Group by protein for efficient PDB loading
    gout["Glutathione_Site"] = gout["Glutathione_Site"].astype(int)
    cyst_by_u = gout.groupby("UniProt_ID")["Glutathione_Site"].apply(list).to_dict()
    ptm_by_u = ptm.groupby("UniProt_ID")

    rows = []
    for uid, cyst_positions in cyst_by_u.items():
        pdb_path = find_alphafold_pdb(uid, pdb_dir)
        if pdb_path is None or not pdb_path.exists():
            continue
        res_atoms = load_residue_atoms(pdb_path)
        # Skip proteins where we cannot map residues
        if not res_atoms:
            continue
        # Build cysteine set: glut + optional controls
        cys_positions = set(int(x) for x in cyst_positions)
        labels_map = {int(x): "glut" for x in cyst_positions}
        if args.include_controls:
            letters = load_residue_letters(pdb_path)
            for pos, aa in letters.items():
                if aa == 'C' and pos not in cys_positions:
                    cys_positions.add(pos)
                    labels_map[pos] = "control"
        df_u = ptm_by_u.get_group(uid) if uid in ptm_by_u.groups else None
        if df_u is None or df_u.empty:
            continue
        for cpos in sorted(cys_positions):
            a_atoms = res_atoms.get(int(cpos))
            if not a_atoms:
                continue
            for _, r in df_u.iterrows():
                ppos = int(r["Position"])
                b_atoms = res_atoms.get(ppos)
                if not b_atoms:
                    continue
                d = min_residue_distance(a_atoms, b_atoms)
                rows.append(
                    {
                        "UniProt_ID": uid,
                        "cysteine_pos": int(cpos),
                        "label": labels_map.get(int(cpos), "glut"),
                        "ptm_type": r["PTM"],
                        "ptm_pos": ppos,
                        "min_distance_A": float(d),
                        "d_le_5A": int(not math.isnan(d) and d <= 5.0),
                        "d_le_8A": int(not math.isnan(d) and d <= 8.0),
                        "d_le_10A": int(not math.isnan(d) and d <= 10.0),
                    }
                )

    out_df = pd.DataFrame(rows)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, index=False)
    print(f"Wrote {len(out_df)} rows -> {args.out}")


if __name__ == "__main__":
    raise SystemExit(main())
