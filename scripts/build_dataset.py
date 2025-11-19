#!/usr/bin/env python3
import argparse
import csv
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd


PTM_FILES = {
    "Acetylation": "Dataset/Acetylation",
    "Methylation": "Dataset/Methylation",
    "Phosphorylation": "Dataset/Phosphorylation",
    "Ubiquitination": "Dataset/Ubiquitination",
    "Succinylation": "Dataset/Succinylation",
    "Malonylation": "Dataset/Malonylation",
    "O_linked_Glycosylation": "Dataset/O-linked Glycosylation",
    "N_linked_Glycosylation": "Dataset/N-linked Glycosylation",
}


# Allowed residue(s) for each PTM based on biology
PTM_ALLOWED_RESIDUE = {
    "Phosphorylation": set("STY"),
    "Ubiquitination": {"K"},
    "Acetylation": {"K"},  # also N-terminus possible, but dataset typically K-centric
    "Methylation": set("KR"),
    "Succinylation": {"K"},
    "Malonylation": {"K"},
    "O_linked_Glycosylation": set("ST"),
    "N_linked_Glycosylation": {"N"},
}


def center_residue(peptide: str) -> str:
    if not peptide:
        return ""
    idx = len(peptide) // 2
    try:
        aa = peptide[idx].upper()
    except Exception:
        return ""
    return aa


def validate_ptm_residue(ptm: str, peptide: str) -> bool:
    aa = center_residue(peptide)
    if aa == "" or aa == "-":
        return False
    allowed = PTM_ALLOWED_RESIDUE[ptm]
    if aa not in allowed:
        return False
    # Optional additional motif check for N-linked: N-X-[S/T] with X != P
    if ptm == "N_linked_Glycosylation":
        i = len(peptide) // 2
        # Need i+1 and i+2 positions
        if i + 2 < len(peptide):
            x = peptide[i + 1].upper()
            z = peptide[i + 2].upper()
            if x == "-" or z == "-":
                return False
            if x == "P":
                return False
            if z not in {"S", "T"}:
                return False
        else:
            return False
    return True


def read_glut_ids(path: Path) -> Tuple[List[str], Dict[str, List[int]]]:
    df = pd.read_csv(path)
    # Expect columns: UniProt_ID, Glutathione_Site
    uniprot_col = None
    for c in df.columns:
        if "uniprot" in c.lower() or "accession" in c.lower():
            uniprot_col = c
            break
    if uniprot_col is None:
        raise ValueError("Could not find UniProt ID column in glutathionylation file")
    site_col = None
    for c in df.columns:
        if "site" in c.lower():
            site_col = c
            break
    if site_col is None:
        raise ValueError("Could not find site column in glutathionylation file")

    df[uniprot_col] = df[uniprot_col].astype(str).str.strip()
    df[site_col] = pd.to_numeric(df[site_col], errors="coerce").astype("Int64")
    df = df.dropna(subset=[uniprot_col, site_col])

    ids = sorted(df[uniprot_col].unique().tolist())
    gl_map: Dict[str, List[int]] = defaultdict(list)
    for uid, grp in df.groupby(uniprot_col):
        gl_map[uid] = sorted(int(x) for x in grp[site_col].dropna().astype(int).tolist())
    return ids, gl_map


def iter_ptm_rows(path: Path, target_ids: set) -> Iterable[Tuple[str, int, str, str]]:
    """
    Yield tuples: (UniProt_ID, position, peptide, ptm_name_from_file)
    File columns (TSV, no header):
    0: Entry name, 1: UniProt, 2: Position, 3: PTM type, 4: Evidence, 5: Peptide (flanking sequence)
    """
    col_names = ["entry", "uniprot", "position", "ptm", "evidence", "peptide"]
    # Use chunked reading for large files
    for chunk in pd.read_csv(path, sep="\t", header=None, names=col_names, chunksize=200000, dtype={2: str}, na_filter=False):
        # Filter to target IDs early
        chunk["uniprot"] = chunk["uniprot"].astype(str).str.strip()
        sub = chunk[chunk["uniprot"].isin(target_ids)].copy()
        if sub.empty:
            continue
        # Coerce position to int where possible
        sub["position"] = pd.to_numeric(sub["position"], errors="coerce").astype("Int64")
        sub = sub.dropna(subset=["position"])
        for _, row in sub.iterrows():
            yield row["uniprot"], int(row["position"]), str(row["peptide"]).strip(), str(row["ptm"]).strip()


def build_dataset(out_path: Path) -> None:
    # 1) Focus UniProt IDs from glutathionylation file
    glut_path = Path("Dataset/cleaned_glutathionylation_sites.csv")
    ids, gl_map = read_glut_ids(glut_path)
    id_set = set(ids)

    # 2) Collect PTM sites for those UniProt IDs
    # Aggregation structure per UniProt
    sites: Dict[str, Dict[str, List[int]]] = defaultdict(lambda: defaultdict(list))
    invalid_sites: Dict[str, Dict[str, List[int]]] = defaultdict(lambda: defaultdict(list))

    for ptm_key, file_path in PTM_FILES.items():
        path = Path(file_path)
        if not path.exists():
            print(f"Warning: PTM file not found: {path}")
            continue
        for uid, pos, pep, ptm_in_file in iter_ptm_rows(path, id_set):
            # Cross-check that file PTM matches mapping (some files embed a specific label)
            # We trust file mapping via ptm_key for allowed residue lookup.
            is_valid = validate_ptm_residue(ptm_key, pep)
            if is_valid:
                sites[uid][ptm_key].append(pos)
            else:
                invalid_sites[uid][ptm_key].append(pos)

    # 3) Write one row per UniProt with semicolon-joined site lists
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Prepare header
    ptm_cols = []
    for ptm_key in PTM_FILES.keys():
        ptm_cols.append(f"{ptm_key}_sites")
        ptm_cols.append(f"{ptm_key}_invalid_sites")

    header = ["UniProt_ID", "Glutathionylation_sites"] + ptm_cols

    with out_path.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)

        for uid in ids:
            row = [uid]
            gl_sites = gl_map.get(uid, [])
            row.append(";".join(str(x) for x in sorted(set(gl_sites))))
            for ptm_key in PTM_FILES.keys():
                valid_list = sorted(set(sites[uid].get(ptm_key, [])))
                invalid_list = sorted(set(invalid_sites[uid].get(ptm_key, [])))
                row.append(";".join(str(x) for x in valid_list))
                row.append(";".join(str(x) for x in invalid_list))
            w.writerow(row)

    print(f"Wrote {out_path}")


def main():
    p = argparse.ArgumentParser(description="Build merged dataset filtered by glutathionylation UniProt IDs with residue validation.")
    p.add_argument("--out", default="Dataset/dataset.csv", help="Output CSV path (default: Dataset/dataset.csv)")
    args = p.parse_args()
    build_dataset(Path(args.out))


if __name__ == "__main__":
    main()

