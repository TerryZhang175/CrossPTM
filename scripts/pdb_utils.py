#!/usr/bin/env python3
from __future__ import annotations

"""
Lightweight PDB helpers for AlphaFold models.

Functions
- load_residue_atoms(pdb_path): {resSeq: [(x,y,z), ...]}
- min_residue_distance(resA, resB): float (angstroms)
"""

from pathlib import Path
from typing import Dict, List, Tuple, Optional
import math


def load_residue_atoms(pdb_path: str | Path) -> Dict[int, List[Tuple[float, float, float]]]:
    res_atoms: Dict[int, List[Tuple[float, float, float]]] = {}
    with open(pdb_path, "r") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            try:
                resseq = int(line[22:26])  # Residue sequence number
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
            except Exception:
                continue
            res_atoms.setdefault(resseq, []).append((x, y, z))
    return res_atoms


def min_residue_distance(a_atoms: List[Tuple[float, float, float]], b_atoms: List[Tuple[float, float, float]]) -> float:
    md = float("inf")
    for (x1, y1, z1) in a_atoms:
        for (x2, y2, z2) in b_atoms:
            dx = x1 - x2
            dy = y1 - y2
            dz = z1 - z2
            d2 = dx * dx + dy * dy + dz * dz
            if d2 < md:
                md = d2
    return math.sqrt(md) if md != float("inf") else float("nan")


_AA3_TO_1 = {
    "ALA": "A", "ARG": "R", "ASN": "N", "ASP": "D",
    "CYS": "C", "GLN": "Q", "GLU": "E", "GLY": "G",
    "HIS": "H", "ILE": "I", "LEU": "L", "LYS": "K",
    "MET": "M", "PHE": "F", "PRO": "P", "SER": "S",
    "THR": "T", "TRP": "W", "TYR": "Y", "VAL": "V",
}


def load_residue_letters(pdb_path: str | Path) -> Dict[int, str]:
    letters: Dict[int, str] = {}
    with open(pdb_path, "r") as fh:
        for line in fh:
            if not line.startswith("ATOM"):
                continue
            try:
                resseq = int(line[22:26])
                res3 = line[17:20].strip().upper()
            except Exception:
                continue
            if resseq not in letters:
                letters[resseq] = _AA3_TO_1.get(res3, "?")
    return letters


def find_alphafold_pdb(uniprot_id: str, pdb_dir: str | Path = "alphafold_structures") -> Optional[Path]:
    """
    Locate an AlphaFold PDB file for the given UniProt ID.
    Accepts any suffix after 'AF-<ID>-F1-' (e.g., model_v4, model_v6, etc.).
    """
    base = Path(pdb_dir)
    pattern = f"AF-{uniprot_id}-F1-*.pdb"
    matches = sorted(base.glob(pattern))
    if matches:
        return matches[0]
    fallback = base / f"AF-{uniprot_id}-F1-model_v4.pdb"
    return fallback if fallback.exists() else None
