#!/usr/bin/env python3
"""
Geometry motifs toolkit (triplet motifs around glutathionylation)

Subcommands:
- prepare-sites: Merge PTMs, restrict to glut proteins, validate residue specificity vs AlphaFold, write normalized CSV.
- detect-motifs: Find 3-node motifs (triangle / wedge) with Cα–Cα cutoff and angle tolerances; CSV only.
- enrich: Residue-type–matched permutation test (Glut fixed) for motif classes; CSV only.

Quickstart:
  python scripts/geometry_motifs.py prepare-sites \
    --glut Dataset/cleaned_glutathionylation_sites.csv \
    --ptm-roots Dataset \
    --af-dir alphafold_structures \
    --out outputs/sites.filtered.csv \
    --drops outputs/sites.dropped.csv
  python scripts/geometry_motifs.py detect-motifs \
    --sites outputs/sites.filtered.csv --out outputs/motifs.csv --cutoff 20 --angle-tol 10
  python scripts/geometry_motifs.py enrich \
    --sites outputs/sites.filtered.csv --out outputs/enrichment.csv --permutations 100
"""

from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
from dataclasses import dataclass
from typing import Dict, Iterable, List, Optional, Set, Tuple

# Path setup for local helpers
from pathlib import Path as _Path
_ROOT = _Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from scripts.pdb_utils import find_alphafold_pdb

# Defaults
DEFAULT_CUTOFF_ANG = 20.0  # Å, Cα–Cα distance cutoff
DEFAULT_ANGLE_TOL = 10.0   # degrees


# --- Utilities ---


AA3_TO_1 = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
}


def infer_uniprot_column(df) -> str:
    cols = list(df.columns)
    cand = [
        c
        for c in cols
        if any(k in c.lower() for k in ["uniprot", "accession", "acc", "protein"])
    ]
    if not cand:
        raise ValueError("Could not find UniProt/accession column")
    return cand[0]


SITE_REGEX = re.compile(r"^([A-Za-z])\s*(\d+)$")


def infer_site_columns(df) -> Tuple[Optional[str], Optional[str]]:
    cols = list(df.columns)
    pos_cols = [
        c
        for c in cols
        if any(k in c.lower() for k in ["position", "pos", "site", "loc", "index"])
    ]
    res_cols = [
        c
        for c in cols
        if any(k in c.lower() for k in ["residue", "aa", "letter"]) and df[c].astype(str).str.len().between(1, 2).any()
    ]
    pos_col = pos_cols[0] if pos_cols else None
    res_col = res_cols[0] if res_cols else None
    return pos_col, res_col


def parse_site_cell(cell: str) -> Tuple[Optional[str], Optional[int]]:
    if cell is None:
        return None, None
    s = str(cell).strip()
    m = SITE_REGEX.match(s)
    if m:
        return m.group(1).upper(), int(m.group(2))
    # Try variants like S-123 or embedded patterns
    m = re.search(r"([A-Za-z])\s*-?\s*(\d+)", s)
    if m:
        return m.group(1).upper(), int(m.group(2))
    # Support number-first pattern like 123K, 48Y
    m = re.search(r"(\d+)\s*-?\s*([A-Za-z])", s)
    if m:
        return m.group(2).upper(), int(m.group(1))
    return None, None


def infer_ptm_type_from_path(path: str) -> str:
    ptm = os.path.basename(os.path.dirname(path)) or "Unknown"
    return ptm


def normalize_ptm_name(name: str) -> str:
    n = name.strip().lower().replace("_", " ")
    mapping = {
        "glutathionylation": "Glutathionylation",
        "s-glutathionylation": "Glutathionylation",
        "phosphorylation": "Phosphorylation",
        "ubiquitination": "Ubiquitination",
        "ubiquitylation": "Ubiquitination",
        "acetylation": "Acetylation",
        "methylation": "Methylation",
        "o-linked glycosylation": "O-linked Glycosylation",
        "o-glycosylation": "O-linked Glycosylation",
        "o-glcnac": "O-linked Glycosylation",
        "n-linked glycosylation": "N-linked Glycosylation",
        "sumoylation": "Sumoylation",
        "nitrosylation": "Nitrosylation",
        "palmitoylation": "Palmitoylation",
    }
    return mapping.get(n, name.strip())


PTM_ALLOWED_RESIDUES: Dict[str, Set[str]] = {
    "Glutathionylation": {"C"},
    "Phosphorylation": {"S", "T", "Y"},
    "Ubiquitination": {"K"},
    "Acetylation": {"K"},
    "Methylation": {"K", "R"},
    "O-linked Glycosylation": {"S", "T"},
    "N-linked Glycosylation": {"N"},
    "Sumoylation": {"K"},
    "Nitrosylation": {"C"},
    "Palmitoylation": {"C"},
}


@dataclass
class AFSequence:
    # pos -> (AA1, pLDDT, x, y, z)
    by_pos: Dict[int, Tuple[str, float, float, float, float]]
    path: str


def load_af_sequence(uniprot: str, af_dir: str) -> Optional[AFSequence]:
    pdb_path = find_alphafold_pdb(uniprot, af_dir)
    if not pdb_path:
        return None
    path = str(pdb_path)
    by_pos: Dict[int, Tuple[str, float, float, float, float]] = {}
    try:
        with open(pdb_path, "r") as fh:
            for line in fh:
                if not line.startswith("ATOM"):
                    continue
                atom = line[12:16].strip()
                if atom != "CA":
                    continue
                resn3 = line[17:20].strip().upper()
                res1 = AA3_TO_1.get(resn3)
                if not res1:
                    continue
                try:
                    resseq = int(line[22:26].strip())
                except ValueError:
                    continue
                try:
                    b = float(line[60:66].strip())  # pLDDT stored in B-factor
                except ValueError:
                    b = float("nan")
                try:
                    x = float(line[30:38].strip())
                    y = float(line[38:46].strip())
                    z = float(line[46:54].strip())
                except ValueError:
                    continue
                if resseq not in by_pos:
                    by_pos[resseq] = (res1, b, x, y, z)
    except Exception as e:
        print(f"Warning: failed to read {path}: {e}", file=sys.stderr)
        return None
    return AFSequence(by_pos=by_pos, path=path)


def collect_glut_uniprots(glut_path: str) -> Set[str]:
    import pandas as pd

    d = pd.read_csv(glut_path)
    uc = infer_uniprot_column(d)
    ids = (
        d[uc]
        .astype(str)
        .str.strip()
        .str.replace("\\.\d+$", "", regex=True)
        .unique()
        .tolist()
    )
    return set(ids)


def iter_ptm_files(ptm_roots: List[str], glut_file: str) -> Iterable[str]:
    # Include data files under given roots (CSV/TSV and extensionless tables), excluding the glut file
    glut_abs = os.path.abspath(glut_file)
    for root in ptm_roots:
        for path in glob.glob(os.path.join(root, "**"), recursive=True):
            if os.path.isdir(path):
                continue
            ap = os.path.abspath(path)
            if ap == glut_abs:
                continue
            base = os.path.basename(path)
            if base.lower().endswith((".csv", ".tsv", ".txt")) or base in {
                "Acetylation",
                "Methylation",
                "Phosphorylation",
                "Ubiquitination",
                "O-linked Glycosylation",
                "N-linked Glycosylation",
                "Succinylation",
                "Malonylation",
            }:
                yield path


def _read_table_any(path: str):
    import pandas as pd
    # Peek to detect delimiter
    sample = ""
    try:
        with open(path, "rb") as fh:
            sample = fh.read(4096).decode(errors="ignore")
    except Exception:
        sample = ""
    has_tab = "\t" in sample
    has_comma = "," in sample
    # Heuristic: if tabs present (our PTM data), prefer TSV headerless
    if has_tab and not has_comma:
        try:
            return pd.read_csv(path, sep="\t", header=None)
        except Exception:
            pass
    # Try CSV (with header)
    try:
        return pd.read_csv(path)
    except Exception:
        pass
    # Try TSV as fallback
    try:
        return pd.read_csv(path, sep="\t")
    except Exception:
        pass
    # Auto-detect as last resort
    return pd.read_csv(path, sep=None, engine="python")


def load_and_normalize_sites(
    glut_path: str,
    ptm_roots: List[str],
    af_dir: str,
    out_path: str,
    drops_path: Optional[str] = None,
    log_non_af_drops: bool = False,
) -> None:
    import pandas as pd

    glut_ids = collect_glut_uniprots(glut_path)

    # Build AF sequences for only glut proteins
    af_cache: Dict[str, AFSequence] = {}
    for uid in sorted(glut_ids):
        af = load_af_sequence(uid, af_dir)
        if af is not None:
            af_cache[uid] = af

    if not af_cache:
        raise SystemExit("No AlphaFold PDBs found for glutathionylated proteins.")

    rows_out = []
    rows_drop = []

    # Always include glutathionylation sites from the glut CSV first
    try:
        d_glu = pd.read_csv(glut_path)
    except Exception as e:
        raise SystemExit(f"Failed to read glut file: {e}")
    uc = infer_uniprot_column(d_glu)
    pos_col, res_col = infer_site_columns(d_glu)
    for i, r in d_glu.iterrows():
        uid = str(r[uc]).strip()
        uid = re.sub(r"\\.\d+$", "", uid)
        if uid not in af_cache:
            if log_non_af_drops:
                rows_drop.append({
                    "reason": "no_af_model",
                    "ptm_type": "Glutathionylation",
                    "uniprot": uid,
                    "position": None,
                    "residue": None,
                    "source_file": os.path.basename(glut_path),
                })
            continue
        pos = None
        res = None
        if pos_col is not None:
            try:
                pos = int(r[pos_col])
            except Exception:
                pos = None
        if res_col is not None and pd.notna(r[res_col]):
            rs = str(r[res_col]).strip().upper()
            res = rs[0] if rs else None
        if pos is None or res is None:
            # try to parse from a site-like column if exists
            for c in d_glu.columns:
                aa, pp = parse_site_cell(r.get(c, None))
                if aa and pp:
                    res, pos = aa, pp
                    break
        if pos is None:
            rows_drop.append({
                "reason": "no_position",
                "ptm_type": "Glutathionylation",
                "uniprot": uid,
                "position": None,
                "residue": res,
                "source_file": os.path.basename(glut_path),
            })
            continue
        af = af_cache[uid]
        aa_plddt = af.by_pos.get(pos)
        if aa_plddt is None:
            rows_drop.append({
                "reason": "pos_missing_in_af",
                "ptm_type": "Glutathionylation",
                "uniprot": uid,
                "position": pos,
                "residue": res,
                "source_file": os.path.basename(glut_path),
            })
            continue
        aa1, plddt, _x, _y, _z = aa_plddt
        if res is None:
            res = aa1
        if res != "C":
            rows_drop.append({
                "reason": "glut_not_on_C",
                "ptm_type": "Glutathionylation",
                "uniprot": uid,
                "position": pos,
                "residue": res,
                "source_file": os.path.basename(glut_path),
            })
            continue
        rows_out.append({
            "uniprot": uid,
            "position": pos,
            "residue": res,
            "ptm_type": "Glutathionylation",
            "plddt": plddt,
            "af_path": af.path,
            "source_file": os.path.basename(glut_path),
        })

    # Other PTMs: only on proteins with glutathionylation and with AF model
    for path in iter_ptm_files(ptm_roots, glut_path):
        try:
            df = _read_table_any(path)
        except Exception as e:
            print(f"Warning: skip {path}: {e}", file=sys.stderr)
            continue
        if df.empty:
            continue
        # Map columns: headered or headerless TSV with 6 cols
        if list(df.columns) == list(range(df.shape[1])):
            # headerless; attempt to map typical 6-col schema
            # 0: entry, 1: uniprot, 2: position, 3: ptm_type, 4: refs, 5: peptide
            uc = 1
            pos_col = 2
            res_col = None
            ptm_col = 3
        else:
            uc = infer_uniprot_column(df)
            pos_col, res_col = infer_site_columns(df)
            ptm_col = None
        ptm_raw = infer_ptm_type_from_path(path)
        ptm = normalize_ptm_name(ptm_raw)
        allowed = PTM_ALLOWED_RESIDUES.get(ptm)

        for i, r in df.iterrows():
            uid = str(r[uc]).strip()
            uid = re.sub(r"\\.\d+$", "", uid)
            if uid not in glut_ids:
                if log_non_af_drops:
                    rows_drop.append({
                        "reason": "protein_not_glut",
                        "ptm_type": ptm,
                        "uniprot": uid,
                        "position": None,
                        "residue": None,
                        "source_file": os.path.basename(path),
                    })
                continue
            af = af_cache.get(uid)
            if af is None:
                if log_non_af_drops:
                    rows_drop.append({
                        "reason": "no_af_model",
                        "ptm_type": ptm,
                        "uniprot": uid,
                        "position": None,
                        "residue": None,
                        "source_file": os.path.basename(path),
                    })
                continue
            # determine position and residue
            pos = None
            res = None
            if pos_col is not None and r.get(pos_col) is not None:
                try:
                    pos = int(r[pos_col])
                except Exception:
                    pos = None
            if res_col is not None and r.get(res_col) is not None:
                try:
                    res = str(r[res_col]).strip().upper()[0]
                except Exception:
                    res = None
            if pos is None or res is None:
                # attempt to parse from any site-like column
                for c in df.columns:
                    aa, pp = parse_site_cell(r.get(c, None))
                    if aa and pp:
                        res, pos = aa, pp
                        break
            # If the table encodes PTM type in-column, override inferred type
            if ptm_col is not None and r.get(ptm_col) is not None:
                ptm_in = normalize_ptm_name(str(r[ptm_col]))
                if ptm_in:
                    ptm = ptm_in
                    allowed = PTM_ALLOWED_RESIDUES.get(ptm)
            if pos is None:
                rows_drop.append({
                    "reason": "no_position",
                    "ptm_type": ptm,
                    "uniprot": uid,
                    "position": None,
                    "residue": res,
                    "source_file": os.path.basename(path),
                })
                continue
            aa_plddt = af.by_pos.get(pos)
            if aa_plddt is None:
                rows_drop.append({
                    "reason": "pos_missing_in_af",
                    "ptm_type": ptm,
                    "uniprot": uid,
                    "position": pos,
                    "residue": res,
                    "source_file": os.path.basename(path),
                })
                continue
            aa1, plddt, _x, _y, _z = aa_plddt
            if res is None:
                res = aa1
            # Enforce residue-specificity if we know allowed residues
            if allowed is not None and res not in allowed:
                rows_drop.append({
                    "reason": "bad_residue",
                    "ptm_type": ptm,
                    "uniprot": uid,
                    "position": pos,
                    "residue": res,
                    "source_file": os.path.basename(path),
                })
                continue
            rows_out.append({
                "uniprot": uid,
                "position": pos,
                "residue": res,
                "ptm_type": ptm,
                "plddt": plddt,
                "af_path": af.path,
                "source_file": os.path.basename(path),
            })

    # Deduplicate
    seen = set()
    deduped = []
    for r in rows_out:
        key = (r["uniprot"], r["position"], r["ptm_type"])
        if key in seen:
            continue
        seen.add(key)
        deduped.append(r)

    # Ensure output directory exists
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=[
                "uniprot",
                "position",
                "residue",
                "ptm_type",
                "plddt",
                "af_path",
                "source_file",
            ],
        )
        writer.writeheader()
        writer.writerows(deduped)

    if drops_path:
        os.makedirs(os.path.dirname(drops_path), exist_ok=True)
        with open(drops_path, "w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["reason", "ptm_type", "uniprot", "position", "residue", "source_file"],
            )
            writer.writeheader()
            writer.writerows(rows_drop)

    print(f"Wrote {len(deduped)} sites to {out_path}")
    if drops_path:
        print(f"Wrote {len(rows_drop)} drops to {drops_path}")


def cmd_prepare_sites(args: argparse.Namespace) -> None:
    ptm_roots = args.ptm_roots or ["Dataset"]
    load_and_normalize_sites(
        glut_path=args.glut,
        ptm_roots=ptm_roots,
        af_dir=args.af_dir,
        out_path=args.out,
        drops_path=args.drops,
        log_non_af_drops=args.log_non_af_drops,
    )


# --- Motif detection and enrichment ---

from math import acos, degrees


def _dist(a, b):
    dx = a[0] - b[0]
    dy = a[1] - b[1]
    dz = a[2] - b[2]
    return (dx * dx + dy * dy + dz * dz) ** 0.5


def _angle(p, a, b):
    # angle at point a between vectors a->p and a->b in degrees
    import math

    v1 = (p[0] - a[0], p[1] - a[1], p[2] - a[2])
    v2 = (b[0] - a[0], b[1] - a[1], b[2] - a[2])
    dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]
    n1 = (v1[0] ** 2 + v1[1] ** 2 + v1[2] ** 2) ** 0.5
    n2 = (v2[0] ** 2 + v2[1] ** 2 + v2[2] ** 2) ** 0.5
    if n1 == 0 or n2 == 0:
        return float("nan")
    cosang = max(-1.0, min(1.0, dot / (n1 * n2)))
    return degrees(acos(cosang))


def _triangle_geom_class(angles: Tuple[float, float, float], angle_tol: float) -> str:
    a1, a2, a3 = angles
    # Right if any ~90°
    if any(abs(a - 90.0) <= angle_tol for a in angles):
        return "triangle_right"
    # Equilateral if all ~60°
    if all(abs(a - 60.0) <= angle_tol for a in angles):
        return "triangle_equilateral"
    # Isosceles if any two angles approximately equal
    if abs(a1 - a2) <= angle_tol or abs(a2 - a3) <= angle_tol or abs(a1 - a3) <= angle_tol:
        return "triangle_isosceles"
    return "triangle_scalene"


def _wedge_geom_class(angle_center: float, angle_tol: float) -> str:
    if abs(angle_center - 90.0) <= angle_tol:
        return "wedge_right"
    if angle_center < 80.0:
        return "wedge_acute"
    if angle_center > 100.0:
        return "wedge_obtuse"
    return "wedge_other"


def _load_sites_grouped(sites_csv: str) -> Dict[str, List[Dict]]:
    import pandas as pd

    d = pd.read_csv(sites_csv)
    required = {"uniprot", "position", "residue", "ptm_type", "plddt", "af_path"}
    missing = required - set(map(str, d.columns))
    if missing:
        raise SystemExit(f"Missing required columns in sites CSV: {sorted(missing)}")
    d["position"] = d["position"].astype(int)
    groups: Dict[str, List[Dict]] = {}
    for r in d.to_dict(orient="records"):
        groups.setdefault(r["uniprot"], []).append(r)
    return groups


def _ensure_af_cache(uniprots: Iterable[str], af_dir_default: Optional[str], from_sites: Dict[str, List[Dict]]) -> Dict[str, AFSequence]:
    af_cache: Dict[str, AFSequence] = {}
    for uid in uniprots:
        # Prefer AF path from any site row; else fall back to af_dir_default
        any_path = None
        for r in from_sites.get(uid, []):
            p = r.get("af_path")
            if isinstance(p, str) and p:
                any_path = p
                break
        if any_path and os.path.exists(any_path):
            af_dir = os.path.dirname(any_path)
        else:
            af_dir = af_dir_default or "alphafold_structures"
        af = load_af_sequence(uid, af_dir)
        if af is not None:
            af_cache[uid] = af
    if not af_cache:
        raise SystemExit("No AlphaFold models found for provided sites.")
    return af_cache


def detect_motifs(
    sites_csv: str,
    out_csv: str,
    cutoff: float = DEFAULT_CUTOFF_ANG,
    angle_tol: float = DEFAULT_ANGLE_TOL,
) -> int:
    groups = _load_sites_grouped(sites_csv)
    af_cache = _ensure_af_cache(groups.keys(), None, groups)

    headers = [
        "uniprot",
        "pos1",
        "pos2",
        "pos3",
        "aa1",
        "aa2",
        "aa3",
        "ptm1",
        "ptm2",
        "ptm3",
        "topology",
        "geom_class",
        "dist12",
        "dist23",
        "dist31",
        "ang1",
        "ang2",
        "ang3",
        "mean_plddt",
    ]

    os.makedirs(os.path.dirname(out_csv), exist_ok=True)
    out_rows = []

    for uid, rows in groups.items():
        af = af_cache.get(uid)
        if af is None:
            continue
        # Build site list with coordinates
        sites = []
        for r in rows:
            pos = int(r["position"])   # 1-based
            entry = af.by_pos.get(pos)
            if entry is None:
                continue
            aa1, plddt, x, y, z = entry
            sites.append({
                "pos": pos,
                "aa": aa1,
                "ptm": str(r["ptm_type"]),
                "plddt": float(r.get("plddt", float("nan"))),
                "coord": (x, y, z),
            })
        n = len(sites)
        if n < 3:
            continue
        # Precompute adjacency by cutoff
        adj = [[False] * n for _ in range(n)]
        for i in range(n):
            for j in range(i + 1, n):
                d = _dist(sites[i]["coord"], sites[j]["coord"])
                if d <= cutoff:
                    adj[i][j] = adj[j][i] = True
        # Enumerate triples
        for i in range(n - 2):
            for j in range(i + 1, n - 1):
                for k in range(j + 1, n):
                    # require at least one Glutathionylation in the triple
                    if not (sites[i]["ptm"] == "Glutathionylation" or sites[j]["ptm"] == "Glutathionylation" or sites[k]["ptm"] == "Glutathionylation"):
                        continue
                    e_ij = adj[i][j]
                    e_jk = adj[j][k]
                    e_ik = adj[i][k]
                    ecount = int(e_ij) + int(e_jk) + int(e_ik)
                    if ecount < 2:
                        continue  # ignore 0/1-edge shapes
                    # Sort nodes by positions for canonical output
                    ord_indices = sorted([i, j, k], key=lambda t: sites[t]["pos"])  # stable
                    ia, ib, ic = ord_indices
                    a = sites[ia]
                    b = sites[ib]
                    c = sites[ic]
                    d12 = _dist(a["coord"], b["coord"])
                    d23 = _dist(b["coord"], c["coord"])
                    d31 = _dist(c["coord"], a["coord"])
                    ang1 = _angle(b["coord"], a["coord"], c["coord"])  # at a
                    ang2 = _angle(a["coord"], b["coord"], c["coord"])  # at b
                    ang3 = _angle(a["coord"], c["coord"], b["coord"])  # at c
                    mean_plddt = sum(x for x in [a["plddt"], b["plddt"], c["plddt"]] if isinstance(x, float)) / 3.0

                    if e_ij and e_jk and e_ik:
                        topology = "triangle"
                        geom_class = _triangle_geom_class((ang1, ang2, ang3), angle_tol)
                    else:
                        topology = "wedge"
                        # Determine center (shared node of the two edges)
                        # edges based on original i,j,k; recompute on ia,ib,ic
                        e_ab = d12 <= cutoff
                        e_bc = d23 <= cutoff
                        e_ac = d31 <= cutoff
                        if e_ab and e_bc:
                            center_angle = ang2
                        elif e_ab and e_ac:
                            center_angle = ang1
                        else:
                            center_angle = ang3
                        geom_class = _wedge_geom_class(center_angle, angle_tol)

                    out_rows.append({
                        "uniprot": uid,
                        "pos1": a["pos"],
                        "pos2": b["pos"],
                        "pos3": c["pos"],
                        "aa1": a["aa"],
                        "aa2": b["aa"],
                        "aa3": c["aa"],
                        "ptm1": a["ptm"],
                        "ptm2": b["ptm"],
                        "ptm3": c["ptm"],
                        "topology": topology,
                        "geom_class": geom_class,
                        "dist12": round(d12, 3),
                        "dist23": round(d23, 3),
                        "dist31": round(d31, 3),
                        "ang1": round(ang1, 2) if ang1 == ang1 else "",
                        "ang2": round(ang2, 2) if ang2 == ang2 else "",
                        "ang3": round(ang3, 2) if ang3 == ang3 else "",
                        "mean_plddt": round(mean_plddt, 1) if mean_plddt == mean_plddt else "",
                    })

    with open(out_csv, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=headers)
        writer.writeheader()
        writer.writerows(out_rows)
    print(f"Wrote {len(out_rows)} motifs to {out_csv}")
    return len(out_rows)


def _count_motif_classes(motif_rows: List[Dict]) -> Dict[Tuple[str, str, str], int]:
    counts: Dict[Tuple[str, str, str], int] = {}
    for r in motif_rows:
        topo = r["topology"]
        geom = r["geom_class"]
        ptms = sorted([r["ptm1"], r["ptm2"], r["ptm3"]])
        key = (topo, geom, "|".join(ptms))
        counts[key] = counts.get(key, 0) + 1
    return counts


def _bh_fdr(pvals: List[float]) -> List[float]:
    m = len(pvals)
    indexed = sorted(enumerate(pvals), key=lambda x: x[1])
    q = [0.0] * m
    prev = 1.0
    for rank, (idx, p) in enumerate(indexed, start=1):
        qval = min(prev, p * m / rank)
        q[idx] = qval
        prev = qval
    return q


def enrich_motifs(
    sites_csv: str,
    out_csv: str,
    cutoff: float = 20.0,
    angle_tol: float = 10.0,
    permutations: int = 100,
    seed: int = 0,
) -> None:
    import pandas as pd
    import random

    groups = _load_sites_grouped(sites_csv)
    af_cache = _ensure_af_cache(groups.keys(), None, groups)

    # Build observed motif rows via detect logic (internal)
    tmp_out = os.path.join(os.path.dirname(out_csv) or ".", "_motifs.tmp.csv")
    os.makedirs(os.path.dirname(tmp_out), exist_ok=True)
    detect_motifs(sites_csv, tmp_out, cutoff, angle_tol)
    obs = pd.read_csv(tmp_out).to_dict(orient="records")
    try:
        os.remove(tmp_out)
    except Exception:
        pass
    obs_counts = _count_motif_classes(obs)

    # Precompute per-protein coordinate index and residue-type pools
    per_protein_coords: Dict[str, Dict[int, Tuple[float, float, float]]] = {}
    per_protein_residues: Dict[str, Dict[str, List[int]]] = {}
    per_protein_glut_positions: Dict[str, Set[int]] = {}
    per_protein_ptm_positions: Dict[str, Dict[str, List[int]]] = {}

    for uid, rows in groups.items():
        af = af_cache.get(uid)
        if af is None:
            continue
        coords = {}
        residues_by_letter: Dict[str, List[int]] = {}
        for pos, (aa, _b, x, y, z) in af.by_pos.items():
            coords[pos] = (x, y, z)
            residues_by_letter.setdefault(aa, []).append(pos)
        per_protein_coords[uid] = coords
        per_protein_residues[uid] = residues_by_letter
        # observed PTM positions
        ptm_pos: Dict[str, List[int]] = {}
        gluts: Set[int] = set()
        for r in rows:
            t = str(r["ptm_type"])
            p = int(r["position"])
            ptm_pos.setdefault(t, []).append(p)
            if t == "Glutathionylation":
                gluts.add(p)
        per_protein_ptm_positions[uid] = ptm_pos
        per_protein_glut_positions[uid] = gluts

    random.seed(seed)
    perm_counts_list: List[Dict[Tuple[str, str, str], int]] = []

    # Helper: detect motifs given a generated PTM assignment per protein
    def detect_counts_from_assignment(assignments: Dict[str, Dict[str, List[int]]]) -> Dict[Tuple[str, str, str], int]:
        rows_all: List[Dict] = []
        for uid, ptmpos in assignments.items():
            coords = per_protein_coords.get(uid)
            if not coords:
                continue
            # Build site list for this permutation
            sites = []
            for t, poses in ptmpos.items():
                for p in poses:
                    coord = coords.get(p)
                    if coord is None:
                        continue
                    # AA letter from residues map
                    aa = None
                    for aa_letter, plist in per_protein_residues[uid].items():
                        # quick membership check; for speed, better build reverse map
                        # but proteins are small, acceptable
                        if p in plist:
                            aa = aa_letter
                            break
                    if aa is None:
                        continue
                    sites.append({"pos": p, "aa": aa, "ptm": t, "coord": coord, "plddt": 0.0})
            n = len(sites)
            if n < 3:
                continue
            adj = [[False] * n for _ in range(n)]
            for i in range(n):
                for j in range(i + 1, n):
                    if _dist(sites[i]["coord"], sites[j]["coord"]) <= cutoff:
                        adj[i][j] = adj[j][i] = True
            # triples
            for i in range(n - 2):
                for j in range(i + 1, n - 1):
                    for k in range(j + 1, n):
                        if not (sites[i]["ptm"] == "Glutathionylation" or sites[j]["ptm"] == "Glutathionylation" or sites[k]["ptm"] == "Glutathionylation"):
                            continue
                        e_ij = adj[i][j]
                        e_jk = adj[j][k]
                        e_ik = adj[i][k]
                        ecount = int(e_ij) + int(e_jk) + int(e_ik)
                        if ecount < 2:
                            continue
                        # Canonical order by positions
                        ord_indices = sorted([i, j, k], key=lambda t: sites[t]["pos"])  # stable
                        ia, ib, ic = ord_indices
                        a = sites[ia]
                        b = sites[ib]
                        c = sites[ic]
                        d12 = _dist(a["coord"], b["coord"])
                        d23 = _dist(b["coord"], c["coord"])
                        d31 = _dist(c["coord"], a["coord"])
                        ang1 = _angle(b["coord"], a["coord"], c["coord"])  # at a
                        ang2 = _angle(a["coord"], b["coord"], c["coord"])  # at b
                        ang3 = _angle(a["coord"], c["coord"], b["coord"])  # at c
                        if e_ij and e_jk and e_ik:
                            topo = "triangle"
                            geom = _triangle_geom_class((ang1, ang2, ang3), angle_tol)
                        else:
                            topo = "wedge"
                            e_ab = d12 <= cutoff
                            e_bc = d23 <= cutoff
                            e_ac = d31 <= cutoff
                            if e_ab and e_bc:
                                center_angle = ang2
                            elif e_ab and e_ac:
                                center_angle = ang1
                            else:
                                center_angle = ang3
                            geom = _wedge_geom_class(center_angle, angle_tol)
                        key = (topo, geom, "|".join(sorted([a["ptm"], b["ptm"], c["ptm"]])))
                        rows_all.append({"key": key})
        # count
        counts: Dict[Tuple[str, str, str], int] = {}
        for r in rows_all:
            k = r["key"]
            counts[k] = counts.get(k, 0) + 1
        return counts

    # Build observed assignment (for consistency in permutation framework)
    obs_assignments: Dict[str, Dict[str, List[int]]] = {}
    for uid, ptmpos in per_protein_ptm_positions.items():
        obs_assignments[uid] = {t: list(set(ps)) for t, ps in ptmpos.items()}

    # Precompute counts per PTM type (excluding Glut fixed)
    def make_perm_assignment() -> Dict[str, Dict[str, List[int]]]:
        assign: Dict[str, Dict[str, List[int]]] = {}
        for uid in groups.keys():
            residues = per_protein_residues.get(uid)
            if not residues:
                continue
            ptmpos = per_protein_ptm_positions.get(uid, {})
            this: Dict[str, List[int]] = {}
            # Keep glut fixed
            glut_positions = per_protein_glut_positions.get(uid, set())
            if glut_positions:
                this["Glutathionylation"] = sorted(glut_positions)
            # For each other PTM, sample positions from allowed residues
            for ptm, allowed in PTM_ALLOWED_RESIDUES.items():
                if ptm == "Glutathionylation":
                    continue
                obs_positions = ptmpos.get(ptm, [])
                need = len(obs_positions)
                if need == 0:
                    continue
                # candidate pool = union of residues with allowed letters in AF
                pool: List[int] = []
                for aa in allowed:
                    pool.extend(residues.get(aa, []))
                if not pool:
                    # no candidates; allocate none
                    continue
                # sample with replacement to allow co-occupancy by different PTMs
                chosen = [random.choice(pool) for _ in range(need)]
                this[ptm] = chosen
            assign[uid] = this
        return assign

    # Run permutations
    perm_counts_list = []
    for _ in range(permutations):
        perm_assign = make_perm_assignment()
        perm_counts = detect_counts_from_assignment(perm_assign)
        perm_counts_list.append(perm_counts)

    # Aggregate permutation counts per key
    keys = set(obs_counts.keys())
    for pc in perm_counts_list:
        keys.update(pc.keys())
    keys = sorted(keys)

    # Compute stats
    data = []
    for k in keys:
        obs_c = obs_counts.get(k, 0)
        perm_values = [pc.get(k, 0) for pc in perm_counts_list]
        mean_exp = sum(perm_values) / len(perm_values) if perm_values else 0.0
        # empirical one-sided p-value (>=)
        ge = sum(1 for v in perm_values if v >= obs_c)
        pval = (ge + 1) / (len(perm_values) + 1)
        data.append({"key": k, "observed": obs_c, "expected": mean_exp, "p_value": pval})

    # FDR
    pvals = [r["p_value"] for r in data]
    qvals = _bh_fdr(pvals) if pvals else []
    for r, q in zip(data, qvals):
        r["q_value"] = q

    # Write CSV
    os.makedirs(os.path.dirname(out_csv) or ".", exist_ok=True)
    with open(out_csv, "w", newline="") as fh:
        writer = csv.writer(fh)
        writer.writerow(["topology", "geom_class", "ptm_multiset", "observed", "expected", "p_value", "q_value"])
        for r in data:
            topo, geom, ptmset = r["key"]
            writer.writerow([topo, geom, ptmset, r["observed"], f"{r['expected']:.3f}", f"{r['p_value']:.4g}", f"{r.get('q_value', float('nan')):.4g}"])


def cmd_detect_motifs(args: argparse.Namespace) -> None:
    detect_motifs(args.sites, args.out, cutoff=args.cutoff, angle_tol=args.angle_tol)


def cmd_enrich(args: argparse.Namespace) -> None:
    enrich_motifs(
        sites_csv=args.sites,
        out_csv=args.out,
        cutoff=args.cutoff,
        angle_tol=args.angle_tol,
        permutations=args.permutations,
        seed=args.seed,
    )


def build_argparser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Geometry motifs toolkit")
    sub = p.add_subparsers(dest="cmd", required=True)

    p_prep = sub.add_parser(
        "prepare-sites", help="Merge PTMs, filter to glut proteins, validate residues via AlphaFold"
    )
    p_prep.add_argument("--glut", required=True, help="Path to cleaned_glutathionylation_sites.csv")
    p_prep.add_argument(
        "--ptm-roots",
        nargs="*",
        default=["Dataset"],
        help="Root directories to scan for PTM CSV files",
    )
    p_prep.add_argument("--af-dir", default="alphafold_structures", help="AlphaFold PDB directory")
    p_prep.add_argument("--out", required=True, help="Output CSV path for normalized sites")
    p_prep.add_argument("--drops", help="Optional CSV path for dropped records")
    p_prep.add_argument("--log-non-af-drops", action="store_true", help="Log drops for proteins not in Glut set or without AF model")
    p_prep.set_defaults(func=cmd_prepare_sites)

    p_det = sub.add_parser("detect-motifs", help="Detect 3-node motifs (triangle/wedge) and write motifs CSV")
    p_det.add_argument("--sites", required=True, help="Normalized sites CSV from prepare-sites")
    p_det.add_argument("--out", required=True, help="Output motifs CSV path")
    p_det.add_argument("--cutoff", type=float, default=DEFAULT_CUTOFF_ANG, help="Cα–Cα distance cutoff (Å)")
    p_det.add_argument("--angle-tol", type=float, default=DEFAULT_ANGLE_TOL, help="Angle tolerance (degrees)")
    p_det.set_defaults(func=cmd_detect_motifs)

    p_enr = sub.add_parser("enrich", help="Permutation enrichment of motif classes (CSV only)")
    p_enr.add_argument("--sites", required=True, help="Normalized sites CSV from prepare-sites")
    p_enr.add_argument("--out", required=True, help="Output enrichment CSV path")
    p_enr.add_argument("--cutoff", type=float, default=DEFAULT_CUTOFF_ANG, help="Cα–Cα distance cutoff (Å)")
    p_enr.add_argument("--angle-tol", type=float, default=DEFAULT_ANGLE_TOL, help="Angle tolerance (degrees)")
    p_enr.add_argument("--permutations", type=int, default=100, help="Number of residue-type–matched shuffles")
    p_enr.add_argument("--seed", type=int, default=0, help="Random seed")
    p_enr.set_defaults(func=cmd_enrich)

    return p


def main(argv: Optional[List[str]] = None) -> int:
    ap = build_argparser()
    args = ap.parse_args(argv)
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
