# AlphaFold Coverage Report

## Scope
- Targets: UniProt IDs from `Dataset/cleaned_glutathionylation_sites.csv` (glutathionylated proteins only).
- Structures: AlphaFold models named `AF-<UniProt>-F1-*.pdb` (e.g., `model_v4`) under `alphafold_structures/`.

## Method
- Extracted target IDs from column `UniProt_ID`.
- Enumerated local AlphaFold files and parsed UniProt accessions.
- For IDs missing locally, queried AlphaFold file endpoints (HTTP HEAD) to classify as available vs. no model.
- Downloaded all available missing models via `scripts/fetch_alphafold_structures.py`.

## Results
- Targets: 694 proteins.
- Local models before fetch: 602.
- Missing before fetch: 92.
- Remote availability of the 92: 62 available, 30 no model (404).
- Local models after fetch: 664.
- Remaining missing (no model): 30.

## Artifacts
- Targets: `reports/artifacts/ids.glut.txt`.
- Local models (after fetch): `reports/artifacts/ids.af.txt`.
- Pre‑fetch missing: `reports/artifacts/ids.missing.txt`.
- Available remotely: `reports/artifacts/ids.available.txt`.
- No model on AlphaFold: `reports/artifacts/ids.no_model.txt`.
- Final missing (after fetch): `reports/artifacts/ids.missing_after.txt`.

## Notes
- PDBs are not tracked in Git; only the manifest `alphafold_structures/download_summary.json` is versioned.
- Reproduce downloads: `python scripts/fetch_alphafold_structures.py --manifest alphafold_structures/download_summary.json --out alphafold_structures`.
- Quote paths with spaces (e.g., `Dataset/O-linked Glycosylation`).
