#!/usr/bin/env python3
"""
Compute residue SASA for selected residue types, preferring FreeSASA if
available and falling back to a built-in Shrake–Rupley implementation.

Usage:
  python scripts/compute_sasa.py --ids-file reports/artifacts/ids.af.txt \
    --pdb-dir alphafold_structures --residue-list C \
    --out reports/sasa_cys.csv

Notes
- If Python freesasa is available, it is used for accuracy and speed.
- Otherwise, a lightweight Shrake–Rupley method is used (sufficient for
  per-residue SASA, especially when limited to a small residue set like C).
"""
from __future__ import annotations

import argparse
from pathlib import Path
import csv
import shutil
import sys
import math
from typing import Dict, List, Tuple

# Local imports for PDB parsing helpers
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from scripts.pdb_utils import load_residue_atoms, load_residue_letters, find_alphafold_pdb

# Local 3-letter to 1-letter mapping (matches scripts/pdb_utils.py)
_AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def have_python_freesasa() -> bool:
    try:
        import freesasa  # noqa: F401
        return True
    except Exception:
        return False


def compute_sasa_python(pdb_path: Path, residue_filter: set[str] | None):
    import freesasa

    struct = freesasa.Structure(str(pdb_path))
    result = freesasa.calc(struct)
    # Build per-residue SASA by iterating atoms
    res_sasa = {}
    for i in range(struct.nAtoms()):
        res_name = struct.residueName(i).strip().upper()
        res_num = struct.residueNumber(i)
        chain = struct.chainLabel(i)
        # Map 3-letter residue names to 1-letter codes
        aa_letter = _AA3_TO_1.get(res_name, "?")
        if residue_filter and aa_letter not in residue_filter:
            continue
        key = (chain, res_num)
        res_sasa[key] = res_sasa.get(key, 0.0) + result.atomArea(i)
    return res_sasa


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ids-file", required=True, help="File with UniProt IDs")
    ap.add_argument("--pdb-dir", default="alphafold_structures")
    ap.add_argument("--residue-list", default="C", help="Residues to include, e.g. C or STY")
    ap.add_argument("--out", default="reports/sasa.csv")
    ap.add_argument("--positions-csv", default="", help="Optional CSV with columns UniProt_ID,cysteine_pos to restrict residues")
    args = ap.parse_args()

    residue_filter = set([x for x in args.residue_list.strip().upper() if x.isalpha()])
    pdb_dir = Path(args.pdb_dir)
    ids = [l.strip() for l in Path(args.ids_file).read_text().splitlines() if l.strip()]

    use_python = have_python_freesasa()
    have_cli = shutil.which("freesasa") is not None
    # Built-in fallback when neither Python nor CLI is available
    use_builtin = not use_python and not have_cli

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    restrict = None
    if args.positions_csv:
        import pandas as pd
        dfp = pd.read_csv(args.positions_csv)
        # Accept column name variants
        pos_col = 'cysteine_pos' if 'cysteine_pos' in dfp.columns else 'Position' if 'Position' in dfp.columns else None
        id_col = 'UniProt_ID'
        if pos_col is None or id_col not in dfp.columns:
            raise SystemExit("positions CSV must have UniProt_ID and cysteine_pos/Position")
        restrict = {}
        for _, r in dfp.iterrows():
            restrict.setdefault(str(r[id_col]).strip(), set()).add(int(r[pos_col]))

    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["UniProt_ID", "chain", "resnum", "SASA", "Residue"])
        for uid in ids:
            pdb_path = find_alphafold_pdb(uid, pdb_dir)
            if pdb_path is None or not pdb_path.exists():
                continue
            if use_python:
                res_sasa = compute_sasa_python(pdb_path, residue_filter)
            elif use_builtin:
                res_sasa = compute_sasa_builtin(pdb_path, residue_filter)
            else:
                print(f"Skipping {uid}: freesasa CLI parsing not implemented", file=sys.stderr)
                continue
            # Emit rows
            for (chain, resnum), sasa in sorted(res_sasa.items(), key=lambda x: x[0][1]):
                if restrict and (uid not in restrict or int(resnum) not in restrict[uid]):
                    continue
                # Try to resolve residue letter from PDB
                letters = load_residue_letters(pdb_path)
                aa = letters.get(int(resnum), "?")
                w.writerow([uid, chain, resnum, f"{sasa:.3f}", aa])
    print(f"Wrote SASA -> {out_path}")
    return 0


# ---------------- Built-in Shrake–Rupley ----------------

_VDW = {
    'H': 1.20, 'C': 1.70, 'N': 1.55, 'O': 1.52, 'F': 1.47, 'P': 1.80,
    'S': 1.80, 'CL': 1.75, 'BR': 1.85, 'SE': 1.90, 'ZN': 1.39,
}
_PROBE = 1.4


def _element_from_line(line: str) -> str:
    e = line[76:78].strip().upper()
    if e:
        return e
    name = line[12:16].strip().upper()
    # Heuristic: first letter (or two for halogens)
    if name.startswith(('CL', 'BR')):
        return name[:2]
    return name[:1]


def _load_all_atoms(pdb_path: Path) -> List[Tuple[float, float, float, float]]:
    atoms = []  # (x,y,z,r)
    with open(pdb_path, 'r') as fh:
        for line in fh:
            if not line.startswith('ATOM'):
                continue
            try:
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except Exception:
                continue
            ele = _element_from_line(line)
            r = _VDW.get(ele, 1.7)
            atoms.append((x, y, z, r))
    return atoms


def _fibonacci_sphere(n: int) -> List[Tuple[float, float, float]]:
    pts = []
    phi = (1 + 5 ** 0.5) / 2
    for i in range(n):
        z = 1 - (2 * i + 1) / n
        r = (1 - z * z) ** 0.5
        theta = 2 * math.pi * i / phi
        x = r * math.cos(theta)
        y = r * math.sin(theta)
        pts.append((x, y, z))
    return pts


def compute_sasa_builtin(pdb_path: Path, residue_filter: set[str] | None, n_points: int = 64):
    # Build all atoms and per-residue atoms
    all_atoms = _load_all_atoms(pdb_path)
    res_atoms = load_residue_atoms(pdb_path)
    letters = load_residue_letters(pdb_path)
    surf_pts = _fibonacci_sphere(n_points)

    # Build a flat list mapping residue->atom indices in all_atoms
    res_atom_indices = {}
    idx = 0
    coords = []
    radii = []
    with open(pdb_path, 'r') as fh:
        for line in fh:
            if not line.startswith('ATOM'):
                continue
            try:
                resseq = int(line[22:26])
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            except Exception:
                continue
            res_atom_indices.setdefault(resseq, []).append(idx)
            coords.append((x, y, z))
            radii.append(_VDW.get(_element_from_line(line), 1.7))
            idx += 1

    res_sasa = {}
    for resnum, atoms in res_atoms.items():
        aa = letters.get(resnum, '?')
        if residue_filter and aa not in residue_filter:
            continue
        # Gather atom indices corresponding to this residue
        idxs = res_atom_indices.get(resnum, [])
        area = 0.0
        # Build residue-level neighbor atom list by scanning all atoms once
        neighbor_idx: List[int] = []
        # Residue atom centers for quick screening
        res_centers = atoms
        for j in range(len(coords)):
            xj, yj, zj = coords[j]
            rj = radii[j]
            # If close to any residue atom center within a liberal cutoff, include
            for (x0, y0, z0) in res_centers:
                dx = x0 - xj; dy = y0 - yj; dz = z0 - zj
                if dx*dx + dy*dy + dz*dz <= (rj + 10.0) ** 2:
                    neighbor_idx.append(j)
                    break
        neighbor_set = set(neighbor_idx)
        for k, (x, y, z) in enumerate(atoms):
            # Map this residue-atom to global atom index approximately by order
            if k >= len(idxs):
                continue
            i = idxs[k]
            ri = radii[i]
            R = ri + _PROBE
            exposed = 0
            # candidate neighbor atoms for this residue
            nbs = neighbor_set
            for (ux, uy, uz) in surf_pts:
                sx = x + R * ux; sy = y + R * uy; sz = z + R * uz
                occluded = False
                for j in nbs:
                    xj, yj, zj = coords[j]
                    rj = radii[j] + _PROBE
                    dx = sx - xj; dy = sy - yj; dz = sz - zj
                    if dx*dx + dy*dy + dz*dz < rj * rj:
                        occluded = True
                        break
                if not occluded:
                    exposed += 1
            frac = exposed / float(n_points)
            area += 4.0 * math.pi * R * R * frac
        res_sasa[(None, resnum)] = area
    return res_sasa

if __name__ == "__main__":
    raise SystemExit(main())
