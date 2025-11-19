# Repository Guidelines

## Project Structure & Module Organization
- `Dataset/`: PTM data grouped by modification (`Phosphorylation/`, `Ubiquitination/`, `"O-linked Glycosylation"/`). Example: `Dataset/cleaned_glutathionylation_sites.csv`.
- `alphafold_structures/`: AlphaFold PDBs (not tracked). Keep only `download_summary.json` in Git.
- `scripts/`: utilities and one-offs (e.g., `fetch_alphafold_structures.py`).
- `tests/`: `pytest` suites; fixtures under `tests/fixtures/`.

## Build, Test, and Development Commands
- Fetch PDBs: `python scripts/fetch_alphafold_structures.py --manifest alphafold_structures/download_summary.json --out alphafold_structures`.
- Count PTM files: `rg --files "Dataset/Phosphorylation" | wc -l`.
- CSV spot‑check: `python -c "import pandas as p; d=p.read_csv('Dataset/cleaned_glutathionylation_sites.csv'); print(d.head()); d.info()"`.
- Format/lint: `black .`; `ruff .` (if configured).
- Run tests: `pytest -q`.

## Coding Style & Naming Conventions
- Python: 4‑space indent; functions/variables `snake_case`; classes `UpperCamelCase`.
- Quote paths with spaces (e.g., `"O-linked Glycosylation"`).
- AlphaFold filenames follow `AF-<UniProt>-F1-*.pdb` (e.g., `model_v4`, `model_v6`). Skip existing files to avoid duplicates.

## Testing Guidelines
- Framework: `pytest`.
- Naming: `tests/test_*.py`; small fixtures in `tests/fixtures/`.
- Cover: parsing/schema checks; file discovery; SASA via FreeSASA; all‑atom min distances; PTM enrichment near glutathionylated cysteines.

## Commit & Pull Request Guidelines
- Commits: imperative, scoped (e.g., `Add PTM distance calculator`, `Filter impossible PTM sites`).
- PRs: include motivation; data deltas (counts before/after filters); methodology (thresholds, radii); sample outputs/plots; link issues; note follow‑ups.

## Data Integrity
- Do not commit PDBs or large binaries. `.gitignore` excludes `alphafold_structures/*` except `download_summary.json`.
- Use consistent naming and reproducible commands; prefer CSV over binaries.

## AlphaFold Coverage Audit (restricted to glutathionylation IDs)
- Extract target IDs: `python -c "import pandas as p; d=p.read_csv('Dataset/cleaned_glutathionylation_sites.csv'); c=[x for x in d.columns if 'uniprot' in x.lower() or 'accession' in x.lower()][0]; print('\n'.join(sorted(set(d[c].astype(str).str.strip()))))" > ids.glut.txt`.
- List local IDs: `ls alphafold_structures/AF-*-F1-*.pdb 2>/dev/null | sed -E 's#.*/AF-([A-Z0-9]+)-F1-.*\.pdb#\1#' | sort -u > ids.af.txt`.
- Find/fetch missing: `comm -23 ids.glut.txt ids.af.txt > ids.missing.txt`; optional `curl` sanity-checks can target any suffix (e.g., `https://alphafold.ebi.ac.uk/files/AF-{}-F1-model_v4.pdb` or `...-F1-v6.pdb`); `python scripts/fetch_alphafold_structures.py --out alphafold_structures --ids $(cat ids.available.txt)` will attempt multiple suffixes automatically.
