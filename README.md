# Crosstalk Analysis Web UI

This repository already contained the glutathionylation crosstalk pipeline (`scripts/ptm_loader.py`, `scripts/compute_distances.py`, `scripts/enrichment.py`, etc.). The new web application lets you drive that pipeline locally with a graphical workflow so that you can load any primary PTM dataset, inspect available secondary PTMs, and run the existing analysis without rewriting code.

## Requirements

Install the Python dependencies once:

```bash
python3 -m pip install -r requirements.txt
```

The requirements file lists the new FastAPI/uvicorn runtime as well as pandas/numpy, which are already used throughout the scripts.

## Start the UI

```bash
uvicorn app.main:app --reload
```

Then open [http://127.0.0.1:8000](http://127.0.0.1:8000) in your browser.

## Workflow

1. **Upload primary PTM sites** – drop a CSV that has UniProt IDs and residue positions. Specify which columns contain those fields (defaults match `Dataset/cleaned_glutathionylation_sites.csv`) and, if needed, provide the “目标氨基酸” filter (e.g., `C` or `S,T,Y`) so any mismatching rows are dropped just like the legacy `ptm_loader` residue checks. After saving to `reports/webui/primary_input.csv`, the app will inspect `alphafold_structures/` and automatically fetch any missing AlphaFold PDBs via `scripts/fetch_alphafold_structures.py` (requires network access).
2. **Select crosstalk PTMs** – choose from the historical PTM files inside `Dataset/` or upload your own PTM CSV/TSV in the “添加新的 Crosstalk PTM” form. Custom datasets are stored under `reports/webui/custom_ptms/` and appended to the loader manifest (with your specified allowed residues for filtering).
3. **Run analysis** – click “运行分析”. The backend sequentially executes the existing scripts (no logic changes) and still writes to the canonical outputs:
   - `reports/artifacts/ptm_sites.filtered.csv`
   - `reports/distances.csv`
   - `reports/enrichment_summary.md` and `reports/cysteine_summary.csv`
4. **Inspect results** – the UI summarizes pair counts, distances, and filtered enrichment text. Download links point to the files above (served read-only via `/reports/...`).
5. **Refresh without rerun** – use “仅刷新结果” to re-apply PTM filters client-side without rerunning the heavy computations.

## Notes

- The backend reuses the CLI scripts through subprocess calls so the scientific logic remains untouched.
- When new UniProt IDs appear in the primary dataset, the service auto-checks `alphafold_structures/` and runs `scripts/fetch_alphafold_structures.py` with `--ids ...` to download any missing models (chunked). The fetcher automatically loops through suffixes `model_v12`, `v12`, …, `model_v1`, `v1`, `model`, `v` so new AlphaFold releases are picked up without code changes; override the order via `--suffixes` or the `ALPHAFOLD_SUFFIXES` environment variable. Missing IDs are automatically re-queued (up to 5 attempts) so once connectivity is restored the downloads resume without re-uploading. If the environment cannot reach AlphaFold, the UI will report which specific IDs failed and why so you can fetch them manually later.
- `scripts/ptm_loader.py` now skips residue inference when a PTM file has no overlapping UniProt IDs, preventing a pandas assignment error that surfaced with the filtered Phosphorylation dataset.
- Uploaded primary PTM data stays under `reports/webui/`; delete that folder to reset the UI state.
