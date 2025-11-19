#!/usr/bin/env python3
import argparse
import csv
import math
from collections import Counter
from itertools import combinations
from pathlib import Path
from typing import Dict, Iterable, List, Sequence, Tuple

import pandas as pd
import sys

from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from scripts.pdb_utils import find_alphafold_pdb


def load_ca_coords(pdb_path: Path) -> Dict[int, Tuple[float, float, float]]:
    coords: Dict[int, Tuple[float, float, float]] = {}
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                resseq = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                continue
            coords[resseq] = (x, y, z)
    return coords


def parse_positions(cell: str) -> List[int]:
    if not isinstance(cell, str) or not cell:
        return []
    out: List[int] = []
    for x in cell.split(";"):
        x = x.strip()
        if not x:
            continue
        try:
            out.append(int(x))
        except ValueError:
            continue
    return out


def angle_between(u: Tuple[float, float, float], v: Tuple[float, float, float]) -> float:
    ux, uy, uz = u
    vx, vy, vz = v
    du = math.sqrt(ux * ux + uy * uy + uz * uz)
    dv = math.sqrt(vx * vx + vy * vy + vz * vz)
    if du == 0 or dv == 0:
        return float("nan")
    dot = ux * vx + uy * vy + uz * vz
    c = max(-1.0, min(1.0, dot / (du * dv)))
    return math.degrees(math.acos(c))


def round_to_tol(x: float, tol: float) -> float:
    return round(x / tol) * tol


def main():
    ap = argparse.ArgumentParser(description="Find common angle motifs among PTM neighbors around glutathionylation sites.")
    ap.add_argument("--input", default="Dataset/dataset_spatial.csv", help="Per-glut spatial CSV from ptm_distance_filter.py")
    ap.add_argument("--out_detail", default="reports/angle_triplets.csv", help="Output CSV of all angle triplets")
    ap.add_argument("--out_summary", default="reports/angle_triplet_summary.csv", help="Output CSV summarizing popular angle combos")
    ap.add_argument("--tolerance", type=float, default=5.0, help="Angle tolerance for binning (degrees, default 5)")
    ap.add_argument("--min_neighbors", type=int, default=3, help="Minimum neighbors to consider (default 3)")
    args = ap.parse_args()

    df = pd.read_csv(args.input)

    # Collect PTM columns (valid sites only)
    ptm_cols = [
        c for c in df.columns
        if c.endswith("_sites")
        and c != "Glutathionylation_sites"
        and not c.endswith("_invalid_sites")
    ]

    detail_rows: List[List[str]] = []
    counts: Counter[str] = Counter()

    for _, row in df.iterrows():
        if str(row.get("structure_available", "")).lower() != "yes":
            continue
        uid = str(row["UniProt_ID"]).strip()
        gl = int(row["Glutathionylation_site"])
        pdb_path = find_alphafold_pdb(uid, "alphafold_structures")
        if not pdb_path or not pdb_path.exists():
            continue
        try:
            coords = load_ca_coords(pdb_path)
        except Exception:
            continue

        gl_ca = coords.get(gl)
        if gl_ca is None:
            continue

        # Gather unique neighbor positions across PTMs that have coordinates
        neighbors: List[int] = []
        seen = set()
        for c in ptm_cols:
            for pos in parse_positions(row.get(c, "")):
                if pos in seen:
                    continue
                if coords.get(pos) is None:
                    continue
                seen.add(pos)
                neighbors.append(pos)

        if len(neighbors) < args.min_neighbors:
            continue

        # For each combination of three neighbors, compute angle triple at the glut site (angles between vectors to neighbors).
        for (i, j, k) in combinations(neighbors, 3):
            vi = (coords[i][0] - gl_ca[0], coords[i][1] - gl_ca[1], coords[i][2] - gl_ca[2])
            vj = (coords[j][0] - gl_ca[0], coords[j][1] - gl_ca[1], coords[j][2] - gl_ca[2])
            vk = (coords[k][0] - gl_ca[0], coords[k][1] - gl_ca[1], coords[k][2] - gl_ca[2])
            aij = angle_between(vi, vj)
            aik = angle_between(vi, vk)
            ajk = angle_between(vj, vk)
            if any(math.isnan(x) for x in (aij, aik, ajk)):
                continue
            angles = sorted([aij, aik, ajk])
            binned = [round_to_tol(x, args.tolerance) for x in angles]
            combo_key = ";".join(f"{x:.1f}" for x in binned)

            detail_rows.append([
                uid,
                str(gl),
                f"{i};{j};{k}",
                f"{angles[0]:.2f};{angles[1]:.2f};{angles[2]:.2f}",
                combo_key,
                str(len(neighbors)),
            ])
            counts[combo_key] += 1

    # Write details
    out_detail = Path(args.out_detail)
    out_detail.parent.mkdir(parents=True, exist_ok=True)
    with out_detail.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["UniProt_ID", "Glutathionylation_site", "neighbor_triplet", "angles_deg_sorted", "binned_combo", "n_neighbors_total"])
        w.writerows(detail_rows)

    # Write summary
    out_summary = Path(args.out_summary)
    out_summary.parent.mkdir(parents=True, exist_ok=True)
    with out_summary.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["binned_combo", "count"])
        for combo, cnt in counts.most_common():
            w.writerow([combo, cnt])

    print(f"Wrote {out_detail} ({len(detail_rows)} triplets) and {out_summary} ({len(counts)} unique combos)")


if __name__ == "__main__":
    main()
