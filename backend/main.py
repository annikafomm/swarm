# backend/main.py

from __future__ import annotations

# Core imports
from enum import Enum
from pathlib import Path
import json
import os
import shutil
import time
from typing import Optional, Dict, Any

# FastAPI / Starlette
from fastapi import FastAPI, File, UploadFile, Form, HTTPException
from fastapi.middleware.cors import CORSMiddleware
from starlette.responses import RedirectResponse

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Base folder for all uploads (created on startup).
BASE_UPLOAD_DIR = Path.cwd() / "uploads"
BASE_UPLOAD_DIR.mkdir(parents=True, exist_ok=True)

# Allow configuring CORS origins via environment variable (comma-separated).
# Example: ALLOWED_ORIGINS="https://myapp.com,https://staging.myapp.com"
ALLOWED_ORIGINS = [
    o.strip()
    for o in os.getenv("ALLOWED_ORIGINS", "http://localhost:4200,http://127.0.0.1:4200").split(",")
    if o.strip()
]

# Optional maximum file size in megabytes (None disables size check).
# Adjust to your needs; set to e.g. 500 for 500 MB or leave as None.
MAX_FILE_MB: Optional[int] = None

# Enumerations for strict validation of form fields.
class Dataset(str, Enum):
    Visium = "Visium"
    Xenium = "Xenium"

class Method(str, Enum):
    Genie3 = "Genie3"
    Sponge = "Sponge"

# -----------------------------------------------------------------------------
# Application
# -----------------------------------------------------------------------------

app = FastAPI()

# CORS: allow only expected frontend origins (configurable).
app.add_middleware(
    CORSMiddleware,
    allow_origins=ALLOWED_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# -----------------------------------------------------------------------------
# Utility helpers
# -----------------------------------------------------------------------------

def _sanitize_filename(name: str) -> str:
    """
    Make a filename safe by allowing only alphanumeric characters and . _ -.
    Everything else is replaced with an underscore.
    """
    return "".join(c if c.isalnum() or c in "._-" else "_" for c in name)


def _ensure_under_max_size(upload: UploadFile) -> None:
    """
    If MAX_FILE_MB is set, check the upload size and raise HTTP 413 if exceeded.
    Note: This reads the in-memory/temporary file; cost is proportional to file size.
    """
    if MAX_FILE_MB is None:
        return
    # Remember current cursor position, go to end to measure, then restore.
    stream = upload.file
    pos = stream.tell()
    stream.seek(0, 2)  # seek to end
    size_bytes = stream.tell()
    stream.seek(pos, 0)  # restore
    size_mb = size_bytes / (1024 * 1024)
    if size_mb > MAX_FILE_MB:
        raise HTTPException(
            status_code=413,
            detail=f"File '{upload.filename}' exceeds size limit of {MAX_FILE_MB} MB (got ~{size_mb:.1f} MB).",
        )


def save_file(upload: Optional[UploadFile], job_dir: Path) -> Optional[str]:
    """
    Persist an uploaded file inside the given job directory.
    Returns the absolute path as string, or None if no file was provided.

    - Ensures the job directory exists.
    - Sanitizes the original filename.
    - Writes the file stream to disk (binary mode).
    """
    if not upload:
        return None

    _ensure_under_max_size(upload)

    job_dir.mkdir(parents=True, exist_ok=True)
    safe_name = _sanitize_filename(upload.filename or "upload.bin")
    dest = job_dir / safe_name

    with dest.open("wb") as f:
        shutil.copyfileobj(upload.file, f)

    # Return as string for JSON serialization
    return str(dest.resolve())


# -----------------------------------------------------------------------------
# Convenience endpoints
# -----------------------------------------------------------------------------

# Simple health check used by load balancers, monitors, or quick manual checks.
@app.get("/health")
def health():
    return {"ok": True}

# Redirect the root URL to the interactive API docs (Swagger UI).
@app.get("/", include_in_schema=False)
def root():
    return RedirectResponse(url="/docs")

# Small JSON index describing key endpoints (for human orientation).
@app.get("/api", include_in_schema=False)
def api_root():
    return {
        "ok": True,
        "message": "API ready",
        "endpoints": {
            "POST /api/upload": "Upload form-data",
            "GET /health": "Health check",
            "GET /docs": "Interactive documentation (Swagger UI)",
        },
    }

# -----------------------------------------------------------------------------
# Upload endpoint
# -----------------------------------------------------------------------------

# Main endpoint for receiving form fields and files (multipart/form-data).
# It creates a dedicated job directory, saves all provided files there,
# stores a config JSON for reproducibility, and returns a summary payload.
@app.post("/api/upload")
async def upload(
    # Required form fields (validated by Enum):
    dataset: Dataset = Form(...),
    method: Method = Form(...),

    # Boolean options (checkboxes on the frontend):
    normalization: bool = Form(...),
    filteringSpatial: bool = Form(...),
    filteringSingleCell: bool = Form(...),
    tangram: bool = Form(...),

    # JSON string with selected scores (parsed below):
    scores: str = Form(...),

    # Files (spatial is required; others are optional):
    spatialFile: UploadFile = File(...),
    singleCellFile: Optional[UploadFile] = File(None),
    precomputedFile: Optional[UploadFile] = File(None),
    spongeNetwork: Optional[UploadFile] = File(None),
    genieFile: Optional[UploadFile] = File(None),
):
    if spatialFile is None:
        raise HTTPException(status_code=400, detail="Spatial file is required")
    
    # 1) Parse the scores JSON
    try:
        scores_obj: Dict[str, Any] = json.loads(scores) if scores else {}
    except Exception:
        raise HTTPException(status_code=400, detail="Field 'scores' must be valid JSON.")

    # 2) Create a unique job directory (timestamp-based)
    job_id = f"job_{int(time.time() * 1000)}"
    job_dir = BASE_UPLOAD_DIR / job_id
    job_dir.mkdir(parents=True, exist_ok=True)

    # 3) Save all files into that directory
    saved_files = {
        "spatialFile": save_file(spatialFile, job_dir),
        "singleCellFile": save_file(singleCellFile, job_dir),
        "precomputedFile": save_file(precomputedFile, job_dir),
        "spongeNetwork": save_file(spongeNetwork, job_dir),
        "genieFile": save_file(genieFile, job_dir),
    }

    # 4) Build response payload
    payload = {
        "ok": True,
        "jobId": job_id,
        "dataset": dataset.value,
        "method": method.value,
        "normalization": normalization,
        "filteringSpatial": filteringSpatial,
        "filteringSingleCell": filteringSingleCell,
        "tangram": tangram,
        "scores": scores_obj,
        "files": saved_files,
    }

    # 5) Persist a copy of the payload next to the uploaded files
    (job_dir / f"{job_id}_config.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")

    # 6) Return a clean JSON response the frontend can consume
    return payload
