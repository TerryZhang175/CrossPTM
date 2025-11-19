#!/usr/bin/env python3
"""
Download AlphaFold PDB models listed in a manifest (download_summary.json).

Usage:
  python scripts/fetch_alphafold_structures.py \
    --manifest alphafold_structures/download_summary.json \
    --out alphafold_structures

Notes:
- Only downloads missing files (skips existing).
- Default base URL matches AlphaFold DB layout. The downloader tries filename suffixes `model_v12`, `v12`, ..., `model_v1`, `v1`, `model`, `v` so it can automatically pick up future releases.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path
from urllib.request import urlopen, Request
from urllib.error import URLError, HTTPError


DEFAULT_BASE = os.environ.get(
    "ALPHAFOLD_BASE_URL", "https://alphafold.ebi.ac.uk/files"
)
def _generate_default_suffixes(max_version: int = 12) -> list[str]:
    suffixes: list[str] = []
    for n in range(max_version, 0, -1):
        suffixes.append(f"model_v{n}")
        suffixes.append(f"v{n}")
    suffixes.append("model")
    suffixes.append("v")
    return suffixes


DEFAULT_SUFFIXES = _generate_default_suffixes()


def build_url(uniprot_id: str, suffix: str, base: str = DEFAULT_BASE) -> str:
    return f"{base}/AF-{uniprot_id}-F1-{suffix}.pdb"


def download_url(url: str, dest: Path, retries: int = 3, backoff: float = 1.5) -> None:
    for attempt in range(retries):
        try:
            req = Request(url, headers={"User-Agent": "crosstalk-analysis-fetch/1.0"})
            with urlopen(req, timeout=60) as r, open(dest, "wb") as f:
                f.write(r.read())
            return
        except (URLError, HTTPError) as e:
            if attempt == retries - 1:
                raise
            time.sleep(backoff ** attempt)


def download_with_suffixes(
    uniprot_id: str,
    out_dir: Path,
    base: str,
    suffixes: list[str],
    retries: int = 3,
    backoff: float = 1.5,
) -> Path:
    last_error: Exception | None = None
    for suffix in suffixes:
        dest = Path(out_dir) / f"AF-{uniprot_id}-F1-{suffix}.pdb"
        url = build_url(uniprot_id, suffix, base)
        try:
            download_url(url, dest, retries=retries, backoff=backoff)
            return dest
        except Exception as exc:
            last_error = exc
            if dest.exists():
                dest.unlink(missing_ok=True)
    if last_error is None:
        raise RuntimeError(f"No suffixes available for {uniprot_id}")
    raise last_error


def parse_suffixes(arg_value: str | None) -> list[str]:
    env_val = os.environ.get("ALPHAFOLD_SUFFIXES")
    raw = arg_value or env_val
    if not raw:
        return DEFAULT_SUFFIXES
    items = [x.strip() for x in raw.split(",")]
    return [x for x in items if x]


def process_id(
    uid: str,
    out_dir: Path,
    base: str,
    suffixes: list[str],
    retries: int = 3,
    backoff: float = 1.5,
) -> str:
    """
    Returns 'downloaded', 'skipped', or 'error' status.
    """
    existing = sorted(out_dir.glob(f"AF-{uid}-F1-*.pdb"))
    if existing:
        return "skipped"
    try:
        download_with_suffixes(
            uid, out_dir, base, suffixes, retries=retries, backoff=backoff
        )
        return "downloaded"
    except Exception as e:
        # We print the error here so it appears in logs immediately
        sys.stderr.write(f"Failed {uid}: {e}\n")
        return "error"


def main(argv=None) -> int:
    p = argparse.ArgumentParser()
    p.add_argument("--manifest", required=True, help="Path to download_summary.json")
    p.add_argument("--out", required=True, help="Output directory for PDB files")
    p.add_argument(
        "--ids",
        nargs="*",
        help="Optional explicit list of UniProt IDs (overrides manifest)",
    )
    p.add_argument("--base", default=DEFAULT_BASE, help="Base URL for downloads")
    p.add_argument(
        "--suffixes",
        default=None,
        help="Comma-separated list of filename suffixes to try",
    )
    p.add_argument(
        "--workers",
        type=int,
        default=4,
        help="Number of concurrent download threads (default: 4)",
    )
    args = p.parse_args(argv)

    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)

    suffixes = parse_suffixes(args.suffixes)

    if args.ids:
        ids = args.ids
    else:
        with open(args.manifest, "r") as f:
            m = json.load(f)
        ids = m.get("successful_ids") or m.get("ids") or []
        if not ids:
            print("No IDs found in manifest.", file=sys.stderr)
            return 1

    downloaded = 0
    skipped = 0
    errors = 0

    import concurrent.futures
    import threading

    lock = threading.Lock()

    # Thread-safe progress printer
    def update_progress(status: str):
        nonlocal downloaded, skipped, errors
        with lock:
            if status == "downloaded":
                downloaded += 1
            elif status == "skipped":
                skipped += 1
            elif status == "error":
                errors += 1
            print(
                f"PROGRESS: downloaded={downloaded} skipped={skipped} errors={errors}",
                flush=True,
            )

    # Use ThreadPoolExecutor for concurrency
    with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as executor:
        # Map each ID to a future
        futures = {
            executor.submit(
                process_id, uid, out_dir, args.base, suffixes
            ): uid
            for uid in ids
        }

        for future in concurrent.futures.as_completed(futures):
            # Retrieve result
            status = future.result()
            update_progress(status)

    print(
        f"Done. downloaded={downloaded} skipped={skipped} errors={errors} -> {out_dir}"
    )
    return 0 if errors == 0 else 2


if __name__ == "__main__":
    raise SystemExit(main())
