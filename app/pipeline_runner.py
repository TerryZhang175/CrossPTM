from __future__ import annotations

import io
import json
import subprocess
import threading
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import Any, Sequence
import re

import pandas as pd

from scripts.ptm_loader import PTM_FILES, parse_site_fields, load_ptm_table
from scripts.pdb_utils import load_residue_letters, find_alphafold_pdb


@dataclass
class PrimaryPTMMeta:
    ptm_name: str
    uniprot_column: str
    site_column: str
    target_residues: str | None = None


@dataclass
class CustomPTMMeta:
    ptm_name: str
    allowed_residues: str | None = None
    filename: str | None = None


class PipelineRunner:
    """Thin wrapper that prepares user inputs and executes the existing scripts."""

    def __init__(self, root: Path | None = None):
        self.root = Path(root or Path(".")).resolve()
        self.reports_dir = self.root / "reports"
        self.work_dir = self.reports_dir / "webui"
        self.work_dir.mkdir(parents=True, exist_ok=True)
        self.primary_csv = self.work_dir / "primary_input.csv"
        self.state_path = self.work_dir / "state.json"
        self.custom_dir = self.work_dir / "custom_ptms"
        self.custom_dir.mkdir(parents=True, exist_ok=True)
        self.custom_manifest = self.custom_dir / "manifest.json"
        if not self.custom_manifest.exists():
            self.custom_manifest.write_text("[]")
        self._fetch_lock = threading.Lock()
        self._fetch_thread: threading.Thread | None = None
        self._fetch_queue: set[str] = set()
        self._fetch_progress: dict[str, Any] = {
            "status": "idle",
            "attempted": 0,
            "downloaded": 0,
            "errors": 0,
            "remaining_missing": 0,
            "failed_ids": [],
            "failed_messages": [],
            "logs": [],
            "queued": 0,
            "total_missing": 0,
        }
        self._fetch_attempt_counts: dict[str, int] = {}
        self._max_fetch_attempts = 5

    # ------------------------------------------------------------------
    # State helpers
    def _load_state(self) -> dict[str, Any]:
        if self.state_path.exists():
            try:
                return json.loads(self.state_path.read_text())
            except json.JSONDecodeError:
                return {}
        return {}

    def _save_state(self, payload: dict[str, Any]) -> None:
        data = self._load_state()
        data.update(payload)
        self.state_path.write_text(json.dumps(data, indent=2))

    def _save_fetch_progress(self) -> None:
        self._save_state({"structure_fetch_progress": self._fetch_progress})

    # ------------------------------------------------------------------
    def list_ptms(self) -> list[dict[str, Any]]:
        """Enumerate PTM sources defined in the historical pipeline."""
        ds_dir = self.root / "Dataset"
        listing = []
        for name in PTM_FILES:
            path = ds_dir / name
            listing.append(
                {
                    "name": name,
                    "path": str(path.relative_to(self.root)) if path.exists() else str(path),
                    "available": path.exists(),
                    "source": "builtin",
                }
            )
        for entry in self._load_custom_ptms():
            path = Path(entry["path"])
            rel = ""
            try:
                rel = str(path.relative_to(self.root))
            except ValueError:
                rel = str(path)
            listing.append(
                {
                    "name": entry["name"],
                    "path": rel,
                    "available": path.exists(),
                    "source": "custom",
                    "allowed_residues": entry.get("allowed_residues"),
                }
            )
        return listing

    # ------------------------------------------------------------------
    def prepare_primary_dataset(
        self, data: bytes, meta: PrimaryPTMMeta
    ) -> dict[str, Any]:
        """Normalize an uploaded PTM CSV so downstream scripts can reuse it."""
        df = pd.read_csv(io.BytesIO(data))
        if meta.uniprot_column not in df.columns:
            raise ValueError(f"Missing UniProt column: {meta.uniprot_column}")
        if meta.site_column not in df.columns:
            raise ValueError(f"Missing site column: {meta.site_column}")

        df = df.copy()
        df[meta.uniprot_column] = df[meta.uniprot_column].astype(str).str.strip()

        residues: list[str | None] = []
        positions: list[int] = []
        uids: list[str] = []
        for _, row in df.iterrows():
            uid = str(row[meta.uniprot_column]).strip()
            site_val = row[meta.site_column]
            aa = None
            pos = self._coerce_position(site_val)
            if pos is None:
                aa, pos = self._parse_site_token(site_val)
            if pos is None:
                continue
            positions.append(int(pos))
            uids.append(uid)
            residues.append(aa)

        if not positions:
            raise ValueError("No valid site positions were found in the uploaded file.")

        out_df = pd.DataFrame({"UniProt_ID": uids, "Glutathione_Site": positions})
        out_df = out_df.dropna(subset=["UniProt_ID"])
        out_df["UniProt_ID"] = out_df["UniProt_ID"].astype(str).str.strip()
        out_df["Glutathione_Site"] = out_df["Glutathione_Site"].astype(int)
        if residues:
            out_df["Residue"] = pd.Series(residues[: len(out_df)], dtype="object")
            out_df["Residue"] = out_df["Residue"].apply(self._normalize_residue_value)
        out_df["Primary_PTM_Name"] = meta.ptm_name.strip() or "primary_ptm"
        out_df = out_df[out_df["UniProt_ID"] != ""]
        if out_df.empty:
            raise ValueError("Primary dataset is empty after filtering invalid rows.")

        allowed = self._parse_residue_filter(meta.target_residues)
        residue_dropped = 0
        allowed_state = sorted(allowed) if allowed else None
        if allowed:
            out_df, residue_dropped = self._apply_residue_filter(out_df, allowed)
            if out_df.empty:
                raise ValueError(
                    "All rows were filtered out because their residues did not match the targeted amino acids."
                )

        out_df.to_csv(self.primary_csv, index=False)
        structure_fetch = self.queue_structure_fetch(set(out_df["UniProt_ID"]))
        self._save_state(
            {
                "primary_csv": str(self.primary_csv.relative_to(self.root)),
                "primary_uploaded": datetime.utcnow().isoformat() + "Z",
                "primary_ptm_name": meta.ptm_name,
                "primary_mapping": {
                    "uniprot_column": meta.uniprot_column,
                    "site_column": meta.site_column,
                    "target_residues": allowed_state,
                },
                "row_count": int(out_df.shape[0]),
                "filtered_by_residue": residue_dropped,
                "structure_fetch": structure_fetch,
            }
        )
        preview = out_df.head(20).to_dict(orient="records")
        return {
            "rows": int(out_df.shape[0]),
            "preview": preview,
            "filtered_by_residue": residue_dropped,
            "target_residues": allowed_state,
            "structure_fetch": structure_fetch,
        }

    # ------------------------------------------------------------------
    def run_pipeline(self, selected_ptms: Sequence[str] | None = None) -> list[str]:
        if not self.primary_csv.exists():
            raise FileNotFoundError("Upload a primary PTM dataset before running the pipeline.")
        ptm_loader_cmd = ["python3", "scripts/ptm_loader.py", "--glut-file", str(self.primary_csv)]
        if selected_ptms:
            for ptm in selected_ptms:
                ptm_loader_cmd.extend(["--ptm", ptm])
        else:
            ptm_loader_cmd.append("--all-ptms")
        custom_manifest = self.custom_manifest if self._has_custom_ptms() else None
        if custom_manifest and custom_manifest.exists():
            ptm_loader_cmd.extend(["--custom-ptms-manifest", str(custom_manifest)])
        commands = [
            ptm_loader_cmd,
            ["python3", "scripts/compute_distances.py", "--glut", str(self.primary_csv)],
            ["python3", "scripts/enrichment.py", "--glut", str(self.primary_csv)],
        ]
        logs: list[str] = []
        for cmd in commands:
            res = subprocess.run(
                cmd,
                cwd=self.root,
                capture_output=True,
                text=True,
            )
            if res.returncode != 0:
                raise RuntimeError(
                    f"Command {' '.join(cmd)} failed:\nSTDOUT:\n{res.stdout}\nSTDERR:\n{res.stderr}"
                )
            if res.stdout:
                logs.append(res.stdout.strip())
        self._save_state({"last_run": datetime.utcnow().isoformat() + "Z"})
        return logs

    # ------------------------------------------------------------------
    def summarize_results(self, ptm_filters: Sequence[str] | None = None) -> dict[str, Any]:
        ptm_filters = [p for p in ptm_filters or [] if p]
        dist_path = self.reports_dir / "distances.csv"
        enrichment_md = self.reports_dir / "enrichment_summary.md"
        ptm_sites = self.reports_dir / "artifacts" / "ptm_sites.filtered.csv"

        summary: dict[str, Any] = {
            "selected_ptms": ptm_filters,
            "distances_path": str(dist_path.relative_to(self.root)) if dist_path.exists() else None,
            "ptm_sites_path": str(ptm_sites.relative_to(self.root)) if ptm_sites.exists() else None,
            "enrichment_path": str(enrichment_md.relative_to(self.root)) if enrichment_md.exists() else None,
        }

        if dist_path.exists():
            df = pd.read_csv(dist_path)
            if ptm_filters:
                df = df[df["ptm_type"].isin(ptm_filters)].copy()
            summary["pair_count"] = int(df.shape[0])
            if not df.empty:
                grouped = (
                    df.groupby("ptm_type")
                    .agg(
                        pair_count=("ptm_type", "size"),
                        mean_distance=("min_distance_A", "mean"),
                        contacts_5A=("d_le_5A", "sum"),
                        contacts_8A=("d_le_8A", "sum"),
                        contacts_10A=("d_le_10A", "sum"),
                    )
                    .reset_index()
                )
                summary["ptm_breakdown"] = grouped.to_dict(orient="records")
            else:
                summary["ptm_breakdown"] = []
        else:
            summary["pair_count"] = 0
            summary["ptm_breakdown"] = []

        if enrichment_md.exists():
            text = enrichment_md.read_text().splitlines()
            summary["enrichment_markdown"] = self._filter_enrichment_lines(text, ptm_filters)
        else:
            summary["enrichment_markdown"] = []

        return summary

    # ------------------------------------------------------------------
    def get_state(self) -> dict[str, Any]:
        state = self._load_state()
        state["primary_ready"] = self.primary_csv.exists()
        if self.primary_csv.exists():
            state["primary_csv"] = str(self.primary_csv.relative_to(self.root))
        state["custom_ptm_count"] = len(self._load_custom_ptms())
        state["structure_fetch_progress"] = self.get_fetch_progress()
        state["missing_structures"] = self.list_missing_structures()
        return state

    def list_missing_structures(self) -> list[str]:
        if not self.primary_csv.exists():
            return []
        try:
            df = pd.read_csv(self.primary_csv, usecols=["UniProt_ID"])
        except Exception:
            return []
        ids = set(df["UniProt_ID"].astype(str).str.strip())
        missing = self._find_missing_structures(ids)
        with self._fetch_lock:
            failed = self._fetch_progress.get("failed_ids", [])
        return sorted(set(missing) | set(failed))

    # ------------------------------------------------------------------
    @staticmethod
    def _filter_enrichment_lines(lines: list[str], ptm_filters: Sequence[str]) -> list[str]:
        if not ptm_filters:
            return lines
        keep = []
        for line in lines:
            if not line.startswith("- "):
                keep.append(line)
                continue
            token = line[2:].split(" ", 1)[0]
            if token == "any" or token in ptm_filters:
                keep.append(line)
        return keep

    @staticmethod
    def _coerce_position(value: Any) -> int | None:
        if pd.isna(value):
            return None
        try:
            # Handles numeric strings and floats
            return int(float(str(value).strip()))
        except (ValueError, TypeError):
            return None

    @staticmethod
    def _parse_site_token(value: Any) -> tuple[str | None, int | None]:
        if value is None or (isinstance(value, float) and pd.isna(value)):
            return None, None
        aa, pos = parse_site_fields(str(value))
        return aa, pos

    @staticmethod
    def _normalize_residue_value(value: Any) -> str | None:
        if value is None:
            return pd.NA
        try:
            if pd.isna(value):
                return pd.NA
        except Exception:
            pass
        s = str(value).strip()
        if not s:
            return pd.NA
        return s.upper()[:1]

    @staticmethod
    def _parse_residue_filter(raw: str | None) -> set[str] | None:
        if raw is None:
            return None
        if isinstance(raw, (list, set, tuple)):
            cleaned = {str(x).strip().upper()[:1] for x in raw if str(x).strip()}
        else:
            tokens = [tok.strip() for tok in str(raw).replace(";", ",").split(",")]
            if len(tokens) == 1 and len(tokens[0]) > 1 and tokens[0].isalnum():
                cleaned = set(tokens[0].upper())
            else:
                cleaned = {tok.upper() for tok in tokens if tok}
        cleaned = {tok for tok in cleaned if len(tok) == 1}
        return cleaned or None

    def _apply_residue_filter(self, df: pd.DataFrame, allowed: set[str]) -> tuple[pd.DataFrame, int]:
        df = df.copy()
        if "Residue" not in df.columns:
            df["Residue"] = pd.NA
        # Normalize to single uppercase letter when available
        df["Residue"] = df["Residue"].apply(self._normalize_residue_value)
        missing_mask = df["Residue"].isna()
        if missing_mask.any():
            self._fill_missing_residues(df, missing_mask)
        keep_mask = df["Residue"].isna() | df["Residue"].isin(allowed)
        dropped = int((~keep_mask).sum())
        filtered = df[keep_mask].copy()
        filtered["Residue"] = filtered["Residue"].apply(self._normalize_residue_value)
        return filtered, dropped

    def _fill_missing_residues(self, df: pd.DataFrame, mask) -> None:
        cache: dict[str, dict[int, str]] = {}
        for idx in df.index[mask]:
            uid = df.at[idx, "UniProt_ID"]
            pos = int(df.at[idx, "Glutathione_Site"])
            aa = self._lookup_residue_letter(uid, pos, cache)
            if aa:
                df.at[idx, "Residue"] = aa

    def _lookup_residue_letter(
        self, uid: str, pos: int, cache: dict[str, dict[int, str]]
    ) -> str | None:
        if uid not in cache:
            pdb_path = find_alphafold_pdb(uid, self.root / "alphafold_structures")
            if pdb_path and pdb_path.exists():
                try:
                    cache[uid] = load_residue_letters(pdb_path)
                except Exception:
                    cache[uid] = {}
            else:
                cache[uid] = {}
        return cache[uid].get(int(pos))

    # ------------------------------------------------------------------
    def register_custom_ptm(self, data: bytes, meta: CustomPTMMeta) -> dict[str, Any]:
        ptm_name = meta.ptm_name.strip()
        if not ptm_name:
            raise ValueError("PTM name cannot be empty.")
        allowed = self._parse_residue_filter(meta.allowed_residues)
        filename = meta.filename or f"{self._slugify(ptm_name)}.csv"
        ext = Path(filename).suffix or ".csv"
        dest = self.custom_dir / f"{self._slugify(ptm_name)}{ext}"
        dest.write_bytes(data)
        # Validate by attempting to parse via existing loader helper
        try:
            load_ptm_table(dest, ptm_name)
        except Exception:
            dest.unlink(missing_ok=True)
            raise
        manifest = self._load_custom_ptms()
        manifest = [entry for entry in manifest if entry.get("name") != ptm_name]
        entry = {
            "name": ptm_name,
            "path": str(dest),
            "allowed_residues": sorted(allowed) if allowed else None,
            "uploaded": datetime.utcnow().isoformat() + "Z",
        }
        manifest.append(entry)
        self._save_custom_ptms(manifest)
        return entry

    def _slugify(self, name: str) -> str:
        slug = re.sub(r"[^A-Za-z0-9._-]+", "_", name.strip())
        return slug or f"ptm_{int(datetime.utcnow().timestamp())}"

    def _load_custom_ptms(self) -> list[dict[str, Any]]:
        if not self.custom_manifest.exists():
            return []
        try:
            return json.loads(self.custom_manifest.read_text())
        except Exception:
            return []

    def _save_custom_ptms(self, entries: list[dict[str, Any]]) -> None:
        self.custom_manifest.write_text(json.dumps(entries, indent=2))

    def _has_custom_ptms(self) -> bool:
        entries = self._load_custom_ptms()
        return bool(entries)

    # ------------------------------------------------------------------
    def queue_structure_fetch(self, ids: set[str]) -> dict[str, Any]:
        ids = {str(i).strip() for i in ids if str(i).strip()}
        if not ids:
            return self.get_fetch_progress()
        with self._fetch_lock:
            self._fetch_queue |= ids
            self._fetch_progress["status"] = "queued"
            self._fetch_progress["queued"] = len(self._fetch_queue)
            self._save_fetch_progress()
            if self._fetch_thread is None or not self._fetch_thread.is_alive():
                self._fetch_thread = threading.Thread(
                    target=self._fetch_worker, daemon=True
                )
                self._fetch_thread.start()
            progress = self._fetch_progress.copy()
            progress["queued"] = len(self._fetch_queue)
            return progress

    def get_fetch_progress(self) -> dict[str, Any]:
        with self._fetch_lock:
            return self._fetch_progress.copy()

    def _fetch_worker(self) -> None:
        while True:
            with self._fetch_lock:
                if not self._fetch_queue:
                    self._fetch_progress = {
                        "status": "idle",
                        "attempted": 0,
                        "downloaded": 0,
                        "errors": 0,
                        "remaining_missing": 0,
                        "failed_ids": [],
                        "failed_messages": [],
                        "logs": [],
                        "queued": 0,
                    }
                    self._save_fetch_progress()
                    self._fetch_thread = None
                    return
                ids = set(self._fetch_queue)
                self._fetch_queue.clear()
                self._fetch_progress.update(
                    {
                        "status": "running",
                        "total": len(ids),
                        "attempted": 0,
                        "downloaded": 0,
                        "errors": 0,
                        "remaining_missing": len(ids),
                        "failed_ids": [],
                        "failed_messages": [],
                        "logs": [],
                        "queued": 0,
                    }
                )
                self._save_fetch_progress()
            try:
                result = self._ensure_structures(ids, progress_cb=self._update_fetch_progress)
            except Exception as exc:  # pragma: no cover
                with self._fetch_lock:
                    self._fetch_progress = {
                        "status": "error",
                        "message": str(exc),
                        "logs": self._fetch_progress.get("logs", []),
                    }
                    self._save_fetch_progress()
                continue
            with self._fetch_lock:
                result["status"] = "completed"
                self._fetch_progress = result
                self._save_fetch_progress()
                retry_candidates = [
                    uid
                    for uid in result.get("remaining_ids", [])
                    if self._fetch_attempt_counts.get(uid, 0) < self._max_fetch_attempts
                ]
            if retry_candidates:
                self.queue_structure_fetch(set(retry_candidates))

    def _update_fetch_progress(self, update: dict[str, Any]) -> None:
        with self._fetch_lock:
            self._fetch_progress.update(update)
            self._save_fetch_progress()

    # ------------------------------------------------------------------
    def _ensure_structures(self, ids: set[str], progress_cb=None) -> dict[str, Any]:
        ids = {str(i).strip() for i in ids if str(i).strip()}
        if not ids:
            return {
                "missing_before": 0,
                "downloaded": 0,
                "errors": 0,
                "remaining_missing": 0,
                "logs": [],
                "failed_ids": [],
                "failed_messages": [],
                "remaining_ids": [],
            }
        missing = self._find_missing_structures(ids)
        if not missing:
            return {
                "missing_before": 0,
                "downloaded": 0,
                "errors": 0,
                "remaining_missing": 0,
                "logs": [],
                "failed_ids": [],
                "failed_messages": [],
                "remaining_ids": [],
            }
        total_missing = len(missing)
        if progress_cb:
            progress_cb(
                {
                    "status": "running",
                    "total": total_missing,
                    "attempted": 0,
                    "downloaded": 0,
                    "errors": 0,
                    "remaining_missing": total_missing,
                    "failed_ids": [],
                    "failed_messages": [],
                }
            )
        fetch_result = self._download_structures(missing, progress_cb=progress_cb)
        remaining = self._find_missing_structures(ids)
        fetch_result["remaining_missing"] = len(remaining)
        fetch_result["missing_before"] = len(missing)
        fetch_result["remaining_ids"] = remaining
        return fetch_result

    def _find_missing_structures(self, ids: set[str]) -> list[str]:
        pdb_dir = self.root / "alphafold_structures"
        pdb_dir.mkdir(parents=True, exist_ok=True)
        missing = []
        for uid in ids:
            pdb_path = find_alphafold_pdb(uid, pdb_dir)
            if not pdb_path or not pdb_path.exists():
                missing.append(uid)
        return missing

    def _download_structures(self, ids: list[str], progress_cb=None) -> dict[str, Any]:
        if not ids:
            return {
                "downloaded": 0,
                "errors": 0,
                "logs": [],
                "attempted": 0,
                "failed_ids": [],
                "failed_messages": [],
            }
        chunk_size = 50
        logs: list[str] = []
        downloaded = 0
        errors = 0
        attempted = 0
        failed_ids: set[str] = set()
        failed_messages: list[str] = []
        total = len(ids)
        manifest = self.root / "alphafold_structures" / "download_summary.json"
        manifest.parent.mkdir(parents=True, exist_ok=True)
        if not manifest.exists():
            manifest.write_text('{"ids": [], "successful_ids": []}\n')
        for chunk in self._chunked(ids, chunk_size):
            cmd = [
                "python3",
                "scripts/fetch_alphafold_structures.py",
                "--manifest",
                str(manifest),
                "--out",
                "alphafold_structures",
                "--ids",
                *chunk,
            ]
            res = subprocess.run(
                cmd,
                cwd=self.root,
                capture_output=True,
                text=True,
            )
            stdout = (res.stdout or "").strip()
            stderr = (res.stderr or "").strip()
            if stdout:
                logs.append(stdout)
            if stderr:
                logs.append(stderr)
                for line in stderr.splitlines():
                    line = line.strip()
                    if not line.startswith("Failed "):
                        continue
                    body = line[len("Failed ") :]
                    if not body:
                        continue
                    uid, _, reason = body.partition(":")
                    uid = uid.strip()
                    if uid:
                        failed_ids.add(uid)
                    failed_messages.append(line)
            match = re.search(r"downloaded=(\d+)\s+skipped=(\d+)\s+errors=(\d+)", stdout)
            if match:
                downloaded += int(match.group(1))
                errors += int(match.group(3))
            attempted += len(chunk)
            for uid in chunk:
                self._fetch_attempt_counts[uid] = self._fetch_attempt_counts.get(uid, 0) + 1
                if uid not in failed_ids and uid in self._fetch_attempt_counts:
                    self._fetch_attempt_counts.pop(uid, None)
            if progress_cb:
                progress_cb(
                    {
                        "attempted": attempted,
                        "downloaded": downloaded,
                        "errors": errors,
                        "failed_ids": sorted(failed_ids),
                        "failed_messages": failed_messages[-10:],
                        "remaining_missing": max(total - attempted, 0),
                    }
                )
            if res.returncode not in (0, 2):
                raise RuntimeError(f"AlphaFold fetch failed: {stdout}\n{stderr}")
        return {
            "downloaded": downloaded,
            "errors": errors,
            "logs": logs,
            "attempted": total,
            "failed_ids": sorted(failed_ids),
            "failed_messages": failed_messages,
        }

    @staticmethod
    def _chunked(seq: list[str], size: int):
        for i in range(0, len(seq), size):
            yield seq[i : i + size]
