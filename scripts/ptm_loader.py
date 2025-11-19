#!/usr/bin/env python3
"""
Collect PTM sites from Dataset/* into a normalized table, filtered to
glutathionylated proteins and positions that are plausible for each PTM.

Output: reports/artifacts/ptm_sites.filtered.csv with columns:
  UniProt_ID, PTM, Residue, Position

Notes
- Attempts to auto-detect delimiter and key columns.
- Recognizes and filters by residue-letter rules per PTM when determinable.
- Restricts to UniProt IDs present in cleaned_glutathionylation_sites.csv
  and with a local AlphaFold structure present.
"""
from __future__ import annotations

import argparse
from pathlib import Path
import pandas as pd
import re
import sys
import json
from pathlib import Path as _Path
# Ensure project root is on sys.path for 'scripts' imports when executed directly
_ROOT = _Path(__file__).resolve().parents[1]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))
from scripts.pdb_utils import load_residue_letters, find_alphafold_pdb


PTM_FILES = [
    "Acetylation",
    "Malonylation",
    "Methylation",
    "N-linked Glycosylation",
    "O-linked Glycosylation",
    "Phosphorylation",
    "Succinylation",
    "Ubiquitination",
]

PTM_ALLOWED_AA = {
    "Acetylation": {"K"},
    "Malonylation": {"K"},
    "Methylation": {"K", "R"},
    "N-linked Glycosylation": {"N"},
    "O-linked Glycosylation": {"S", "T"},
    "Phosphorylation": {"S", "T", "Y"},
    "Succinylation": {"K"},
    "Ubiquitination": {"K"},
}


def infer_uniprot_column(df: pd.DataFrame) -> str | None:
    cols = list(df.columns)
    cand = [
        c
        for c in cols
        if any(k in c.lower() for k in ["uniprot", "accession", "acc", "entry"])
    ]
    return cand[0] if cand else None


def infer_site_columns(df: pd.DataFrame) -> tuple[str | None, str | None]:
    aa_col = None
    pos_col = None
    for c in df.columns:
        lc = c.lower()
        if aa_col is None and any(k in lc for k in ["aa", "res", "amino"]):
            aa_col = c
        if pos_col is None and any(k in lc for k in ["position", "site", "pos"]):
            pos_col = c
    return aa_col, pos_col


def parse_site_fields(site_val: str) -> tuple[str | None, int | None]:
    """Parse strings like 'K123', 'S 45', '123', return (AA?, pos)."""
    if site_val is None or (isinstance(site_val, float) and pd.isna(site_val)):
        return None, None
    s = str(site_val).strip()
    m = re.match(r"^([A-Za-z])\s*-?\s*(\d+)$", s)
    if m:
        return m.group(1).upper(), int(m.group(2))
    m = re.match(r"^(\d+)$", s)
    if m:
        return None, int(m.group(1))
    return None, None


def load_ptm_table(path: Path, ptm: str) -> pd.DataFrame:
    # Try headerless TSV with 6 columns (common format)
    try:
        df = pd.read_csv(path, sep="\t", header=None)
        if df.shape[1] >= 3:
            df = df.rename(columns={1: "UniProt_ID", 2: "Position", 3: "PTM", 5: "Peptide"})
    except Exception:
        df = pd.read_csv(path, sep=None, engine="python")
    # Fast path: headerless inferred columns present
    if {"UniProt_ID", "Position"}.issubset(df.columns):
        ucol = "UniProt_ID"
        aa_col = None
        pos_col = "Position"
        # PTM column may not match filename; override with provided ptm label
        df["PTM"] = ptm
    else:
        ucol = infer_uniprot_column(df)
        aa_col, pos_col = infer_site_columns(df)
    if ucol is None or pos_col is None:
        # Try fallback for common spellings
        for c in df.columns:
            if ucol is None and c.strip() in {"UniProt_ID", "UniProt", "ACC_ID"}:
                ucol = c
            if pos_col is None and c.strip() in {"Position", "Site"}:
                pos_col = c
    if ucol is None or pos_col is None:
        # Minimal fallback: look for any integer-like column for position
        for c in df.columns:
            try:
                if df[c].dropna().astype(int).shape[0] > 0 and pos_col is None:
                    pos_col = c
            except Exception:
                pass
    if ucol is None or pos_col is None:
        raise ValueError(f"Could not infer columns from {path}")

    rows = []
    for _, r in df.iterrows():
        uid = str(r[ucol]).strip()
        aa = None
        pos = None
        # Prefer explicit columns
        if aa_col and (aa_col in df.columns) and pd.notna(r.get(aa_col)):
            try:
                aa = str(r[aa_col]).strip().upper()[:1]
            except Exception:
                pass
        if pos_col and (pos_col in df.columns) and pd.notna(r.get(pos_col)):
            try:
                pos = int(r[pos_col])
            except Exception:
                aa2, pos2 = parse_site_fields(r[pos_col])
                aa = aa or aa2
                pos = pos2
        if pos is None and "Site" in df.columns:
            aa2, pos2 = parse_site_fields(r["Site"])
            aa = aa or aa2
            pos = pos2
        if pos is None:
            continue
        rows.append((uid, ptm, aa, pos))
    out = pd.DataFrame(rows, columns=["UniProt_ID", "PTM", "Residue", "Position"])
    return out


def load_custom_ptms_manifest(manifest_path: str | None):
    if not manifest_path:
        return []
    mp = Path(manifest_path)
    if not mp.exists():
        return []
    try:
        data = json.loads(mp.read_text())
    except Exception as exc:
        print(f"Warning: could not parse custom PTM manifest {manifest_path}: {exc}")
        return []
    entries = []
    for item in data:
        name = item.get("name")
        path = item.get("path")
        if not name or not path:
            continue
        allowed = item.get("allowed_residues")
        if isinstance(allowed, str):
            allowed = [allowed]
        entries.append(
            {
                "name": name,
                "path": path,
                "allowed_residues": allowed,
            }
        )
    return entries


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--dataset", default="Dataset", help="Dataset directory")
    ap.add_argument(
        "--glut-file",
        default="Dataset/cleaned_glutathionylation_sites.csv",
        help="CSV of glutathionylation sites (must include UniProt_ID)",
    )
    ap.add_argument(
        "--out",
        default="reports/artifacts/ptm_sites.filtered.csv",
        help="Output CSV path",
    )
    ap.add_argument(
        "--custom-ptms-manifest",
        default=None,
        help="Optional JSON manifest describing user-added PTM files",
    )
    ap.add_argument(
        "--ptm",
        action="append",
        dest="selected_ptms",
        help="Restrict processing to the specified PTM name (can repeat); use --all-ptms to include built-ins.",
    )
    ap.add_argument(
        "--all-ptms",
        action="store_true",
        help="Include all built-in PTMs (ignored if --ptm is used).",
    )
    args = ap.parse_args()

    ds = Path(args.dataset)
    gout = pd.read_csv(args.glut_file)
    if "UniProt_ID" not in gout.columns:
        raise ValueError("Glutathionylation CSV must have UniProt_ID column")
    target_ids = set(gout["UniProt_ID"].astype(str).str.strip())

    custom_ptms = load_custom_ptms_manifest(args.custom_ptms_manifest)

    # Restrict to proteins with local AlphaFold PDBs
    af_dir = Path("alphafold_structures")
    have_af = set(
        p.name.split("-")[1]
        for p in af_dir.glob("AF-*-F1-*.pdb")
        if p.is_file()
    )
    targets = target_ids & have_af

    selected_set = set(args.selected_ptms or [])

    all_rows = []
    ptm_sources = []
    include_builtin = args.all_ptms or not selected_set
    for fname in PTM_FILES:
        if not include_builtin and fname not in selected_set:
            continue
        if selected_set and fname not in selected_set:
            continue
        ptm_sources.append(
            {
                "name": fname,
                "path": ds / fname,
                "allowed": PTM_ALLOWED_AA.get(fname),
            }
        )
    for entry in custom_ptms:
        if selected_set and entry["name"] not in selected_set:
            continue
        ptm_sources.append(
            {
                "name": entry["name"],
                "path": Path(entry["path"]),
                "allowed": set(entry["allowed_residues"]) if entry.get("allowed_residues") else None,
            }
        )

    for source in ptm_sources:
        fname = source["name"]
        fpath = source["path"]
        allowed_override = source["allowed"]
        if not fpath.exists():
            continue
        try:
            df = load_ptm_table(fpath, fname)
        except Exception as e:
            print(f"Skip {fname}: {e}")
            continue
        # Filter to target proteins
        df = df[df["UniProt_ID"].isin(targets)].copy()
        if df.empty:
            # If no targets remain for this PTM we can skip residue inference entirely.
            all_rows.append(df)
            continue

        # If Residue is unknown, attempt to infer from local PDB residue names
        cache_letters: dict[str, dict[int, str]] = {}
        for uid in df["UniProt_ID"].unique():
            pdb_path = find_alphafold_pdb(uid, "alphafold_structures")
            if pdb_path:
                cache_letters[uid] = load_residue_letters(pdb_path)
        def infer_res(row):
            if pd.notna(row["Residue"]) and str(row["Residue"]).strip():
                return str(row["Residue"]).strip().upper()[:1]
            letters = cache_letters.get(row["UniProt_ID"], {})
            return letters.get(int(row["Position"]), None)
        df["Residue"] = df.apply(infer_res, axis=1)
        # Filter by allowed residues when determinable
        allowed = allowed_override if allowed_override is not None else PTM_ALLOWED_AA.get(fname)
        if allowed is not None:
            df = df[(df["Residue"].isna()) | (df["Residue"].isin(allowed))]
        all_rows.append(df)

    if not all_rows:
        print("No PTM files parsed.")
        return 0
    out_df = pd.concat(all_rows, ignore_index=True)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    out_df.to_csv(args.out, index=False)
    print(f"Wrote {len(out_df)} rows -> {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
