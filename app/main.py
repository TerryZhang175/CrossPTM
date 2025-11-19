from __future__ import annotations

from pathlib import Path
from typing import Annotated

from fastapi import FastAPI, File, Form, HTTPException, Query, UploadFile
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles
from pydantic import BaseModel

from app.pipeline_runner import PipelineRunner, PrimaryPTMMeta, CustomPTMMeta


class RunRequest(BaseModel):
    ptm_types: list[str] | None = None


app = FastAPI(title="Crosstalk Analysis UI", version="0.1.0")
runner = PipelineRunner()

static_dir = Path(__file__).resolve().parent / "static"
if static_dir.exists():
    app.mount("/static", StaticFiles(directory=static_dir), name="static")
reports_dir = Path("reports")
if reports_dir.exists():
    app.mount("/reports", StaticFiles(directory=reports_dir), name="reports")


@app.get("/", response_class=HTMLResponse)
def index():
    index_path = static_dir / "index.html"
    if not index_path.exists():
        raise HTTPException(status_code=500, detail="Missing frontend index.html")
    return index_path.read_text()


@app.get("/api/ptms")
def list_ptms():
    return {"ptms": runner.list_ptms()}


@app.get("/api/state")
def get_state():
    return runner.get_state()


@app.post("/api/primary")
async def upload_primary(
    primary_ptm_name: Annotated[str, Form(..., description="Name of the PTM under study")],
    uniprot_column: Annotated[str, Form()] = "UniProt_ID",
    site_column: Annotated[str, Form()] = "Glutathione_Site",
    target_residues: Annotated[str | None, Form()] = None,
    file: UploadFile = File(...),
):
    data = await file.read()
    meta = PrimaryPTMMeta(
        ptm_name=primary_ptm_name,
        uniprot_column=uniprot_column,
        site_column=site_column,
        target_residues=target_residues,
    )
    try:
        info = runner.prepare_primary_dataset(data, meta)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return {"message": "Primary PTM dataset stored.", "info": info}


@app.post("/api/run")
def run_pipeline(payload: RunRequest):
    try:
        logs = runner.run_pipeline(payload.ptm_types)
    except FileNotFoundError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except RuntimeError as exc:
        raise HTTPException(status_code=500, detail=str(exc)) from exc
    summary = runner.summarize_results(payload.ptm_types or [])
    return {"logs": logs, "summary": summary}


@app.get("/api/results")
def get_results(ptm: list[str] | None = Query(default=None)):
    return runner.summarize_results(ptm)


@app.get("/api/structures/progress")
def structure_progress():
    return {
        "progress": runner.get_fetch_progress(),
        "missing": runner.list_missing_structures(),
        "primary_ready": runner.primary_csv.exists(),
    }


@app.post("/api/custom-ptm")
async def upload_custom_ptm(
    ptm_name: Annotated[str, Form(..., description="Name of the additional PTM")],
    allowed_residues: Annotated[str | None, Form()] = None,
    file: UploadFile = File(...),
):
    data = await file.read()
    meta = CustomPTMMeta(
        ptm_name=ptm_name,
        allowed_residues=allowed_residues,
        filename=file.filename,
    )
    try:
        entry = runner.register_custom_ptm(data, meta)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except Exception as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return {"message": "Custom PTM uploaded.", "ptm": entry}
