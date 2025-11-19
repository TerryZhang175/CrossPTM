#!/usr/bin/env python3
import argparse
import csv
import math
from pathlib import Path
from typing import Dict, List, Tuple

import pandas as pd
import sys

from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from scripts.pdb_utils import find_alphafold_pdb


def load_ca_coords(pdb_path: Path) -> Dict[int, Tuple[float, float, float]]:
    """Parse PDB and return mapping: resSeq -> (x,y,z) for CA atoms.
    Assumes AlphaFold numbering matches UniProt positions.
    """
    coords: Dict[int, Tuple[float, float, float]] = {}
    with open(pdb_path, "r", encoding="utf-8", errors="ignore") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            atom_name = line[12:16]
            if atom_name.strip() != "CA":
                continue
            # residue sequence number (right-aligned)
            try:
                resseq = int(line[22:26].strip())
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                continue
            coords[resseq] = (x, y, z)
    return coords


def dist(a: Tuple[float, float, float], b: Tuple[float, float, float]) -> float:
    return math.dist(a, b)


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
            pass
    return out


def main():
    ap = argparse.ArgumentParser(description="Filter PTM sites by C-alpha distance to glutathionylation sites (< threshold). One row per glut site.")
    ap.add_argument("--input", default="Dataset/dataset.csv", help="Input merged dataset (from build_dataset.py)")
    ap.add_argument("--out", default="Dataset/dataset_spatial.csv", help="Output CSV path")
    ap.add_argument("--threshold", type=float, default=20.0, help="Distance threshold in Angstroms (default: 20.0)")
    args = ap.parse_args()

    df = pd.read_csv(args.input)

    # Identify PTM site columns to filter (valid sites only)
    ptm_cols = [
        c for c in df.columns
        if c.endswith("_sites")
        and not c.endswith("_invalid_sites")
        and c != "Glutathionylation_sites"
    ]

    rows_out: List[List[str]] = []
    header = ["UniProt_ID", "Glutathionylation_site", "structure_available"] + ptm_cols

    for _, row in df.iterrows():
        uid = str(row["UniProt_ID"]).strip()
        gl_sites = parse_positions(row.get("Glutathionylation_sites", ""))
        if not gl_sites:
            continue

        # Load structure once per UniProt
        pdb_path = find_alphafold_pdb(uid, "alphafold_structures")
        coords = None
        structure_available = pdb_path is not None and pdb_path.exists()
        if structure_available:
            try:
                coords = load_ca_coords(pdb_path)
            except Exception:
                coords = None
                structure_available = False

        # Pre-parse PTM positions per PTM
        ptm_positions: Dict[str, List[int]] = {c: parse_positions(row.get(c, "")) for c in ptm_cols}

        for gl in gl_sites:
            # If no structure, emit empty PTM lists
            if not structure_available or coords is None:
                out_row = [uid, str(gl), "no"]
                out_row += ["" for _ in ptm_cols]
                rows_out.append(out_row)
                continue

            gl_ca = coords.get(int(gl))
            if gl_ca is None:
                # No CA for that residue number; emit empty
                out_row = [uid, str(gl), "no_CA"]
                out_row += ["" for _ in ptm_cols]
                rows_out.append(out_row)
                continue

            # For each PTM list, keep those within threshold
            kept: Dict[str, List[int]] = {}
            for c in ptm_cols:
                kept[c] = []
                for pos in ptm_positions.get(c, []):
                    ca = coords.get(int(pos))
                    if ca is None:
                        continue
                    d = dist(gl_ca, ca)
                    if d < args.threshold:
                        kept[c].append(int(pos))

            out_row = [uid, str(gl), "yes"]
            for c in ptm_cols:
                out_row.append(";".join(str(x) for x in sorted(set(kept[c]))))
            rows_out.append(out_row)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        w.writerows(rows_out)

    print(f"Wrote {out_path} with {len(rows_out)} rows")


if __name__ == "__main__":
    main()
