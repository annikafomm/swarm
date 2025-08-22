# backend/main.py

from __future__ import annotations

# ---------------------------------
# Standard Library
# ---------------------------------
import json
import os
import shutil
import time
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Optional
from uuid import UUID, uuid4

import pandas as pd
import scanpy as sc
import uvicorn

# ---------------------------------
# Third-Party (FastAPI / Starlette)
# ---------------------------------
from fastapi import (
    Depends,
    FastAPI,
    File,
    Form,
    HTTPException,
    Response,
    UploadFile,
)
from fastapi.middleware.cors import CORSMiddleware
from fastapi_sessions.backends.implementations import InMemoryBackend
from fastapi_sessions.frontends.implementations import (
    CookieParameters,
    SessionCookie,
)
from fastapi_sessions.session_verifier import SessionVerifier
from pydantic import BaseModel as PydanticBaseModel
from pydantic import EmailStr
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
    for o in os.getenv(
        "ALLOWED_ORIGINS", "http://localhost:4200,http://127.0.0.1:4200"
    ).split(",")
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


class BaseModel(PydanticBaseModel):
    class Config:
        arbitrary_types_allowed = True


class SessionData(BaseModel):
    username: str
    adata: sc.AnnData | None = None


class BasicVerifier(SessionVerifier[UUID, SessionData]):
    def __init__(
        self,
        *,
        identifier: str,
        auto_error: bool,
        backend: InMemoryBackend[UUID, SessionData],
        auth_http_exception: HTTPException,
    ):
        self._identifier = identifier
        self._auto_error = auto_error
        self._backend = backend
        self._auth_http_exception = auth_http_exception

    @property
    def identifier(self):
        return self._identifier

    @property
    def backend(self):
        return self._backend

    @property
    def auto_error(self):
        return self._auto_error

    @property
    def auth_http_exception(self):
        return self._auth_http_exception

    def verify_session(self, model: SessionData) -> bool:
        """If the session exists, it is valid"""
        return True


class AnnDataPath(BaseModel):
    path: str


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

# Uses UUID
cookie_params = CookieParameters(secure=False, httponly=True, samesite="lax")

cookie = SessionCookie(
    cookie_name="cookie",
    identifier="general_verifier",
    auto_error=True,
    secret_key="DONOTUSE",
    cookie_params=cookie_params,
)
backend = InMemoryBackend[UUID, SessionData]()

verifier = BasicVerifier(
    identifier="general_verifier",
    auto_error=True,
    backend=backend,
    auth_http_exception=HTTPException(
        status_code=403, detail="invalid session"
    ),
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
@app.post("/api/upload", dependencies=[Depends(cookie)])
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
    useDefaultLiana: Optional[bool] = Form(None),
    tangramOptions: Optional[str] = Form(None),
    spongeOptions: Optional[str] = Form(None),
    squidpyOptions: Optional[str] = Form(None),
    lianaOptions: Optional[str] = Form(None),
    # Files (spatial is required; others are optional):
    email: EmailStr = Form(...),
    spatialFile: UploadFile = File(...),
    singleCellFile: Optional[UploadFile] = File(None),
    precomputedFile: Optional[UploadFile] = File(None),
    spongeNetworkAnalysis: Optional[UploadFile] = File(None),
    spongeNetworkInteractions: Optional[UploadFile] = File(None),
    genieFile: Optional[UploadFile] = File(None),
    session_data: "SessionData" = Depends(verifier),
):
    if spatialFile is None:
        raise HTTPException(status_code=400, detail="Spatial file is required")

    # 1) Parse the scores JSON
    try:
        scores_obj: Dict[str, Any] = json.loads(scores) if scores else {}
    except Exception:
        raise HTTPException(
            status_code=400, detail="Field 'scores' must be valid JSON."
        )

    # 1b) Option-JSONs sicher parsen
    def _parse_json_field(name: str, val: Optional[str]):
        if val in (None, "", "null"):
            return None
        try:
            return json.loads(val)
        except Exception:
            raise HTTPException(
                status_code=400, detail=f"Field '{name}' must be valid JSON."
            )

    options = {
        "tangram": _parse_json_field("tangramOptions", tangramOptions),
        "sponge": _parse_json_field("spongeOptions", spongeOptions),
        "squidpy": _parse_json_field("squidpyOptions", squidpyOptions),
        "liana": _parse_json_field("lianaOptions", lianaOptions),
    }
    # None-Einträge entfernen
    options = {k: v for k, v in options.items() if v is not None}

    # --- Validation rules ----------------------------------------------------
    # Helper: normalize scores to a set of names (list/obj tolerant)
    def _scores_to_set(obj):
        try:
            if isinstance(obj, dict):
                # allow object like {"LIANA+": true, ...}
                return {k for k, v in obj.items() if v}
            if isinstance(obj, list):
                return set(obj)
            return set()
        except Exception:
            return set()

    scores_set = _scores_to_set(scores_obj)

    sponge_needed = (
        ("SPONGeffects" in scores_set)
        or ("AUCell" in scores_set)
        or (method == Method.Sponge)
    ) and (precomputedFile is None)
    genie_needed = (
        ("VIPER" in scores_set)
        or ("AUCell" in scores_set)
        or ("SPONGeffects" in scores_set)
    ) and (precomputedFile is None)

    # spatial always required (already enforced by FastAPI type); double-check for safety
    if spatialFile is None:
        raise HTTPException(
            status_code=400, detail="Spatial file is required."
        )

    # If tangram true -> singleCell required
    if tangram and (singleCellFile is None):
        raise HTTPException(
            status_code=400,
            detail="Tangram is enabled: singleCell file is required.",
        )

    # SPONGE pair requirement
    if sponge_needed:
        missing = []
        if spongeNetworkAnalysis is None:
            missing.append("spongeNetworkAnalysis")
        if spongeNetworkInteractions is None:
            missing.append("spongeNetworkInteractions")
        if missing:
            raise HTTPException(
                status_code=400,
                detail=f"SPONGE requires both files: {', '.join(missing)} missing.",
            )

    # GENIE/VIPER requirement
    if genie_needed and (genieFile is None):
        raise HTTPException(
            status_code=400,
            detail="VIPER/AUCell/SPONGeffects selected: genieNetwork file is required (unless precomputed is provided).",
        )
    raw_username = session_data.username
    user_safe = _sanitize_filename(raw_username) or "anon"

    job_id = f"job_{int(time.time() * 1000)}"
    job_dir = BASE_UPLOAD_DIR / f"{job_id}_{user_safe}"
    job_dir.mkdir(parents=True, exist_ok=True)

    # 3) Save all files into that directory
    saved_files = {
        "spatialFile": save_file(spatialFile, job_dir),
        "singleCellFile": save_file(singleCellFile, job_dir),
        "precomputedFile": save_file(precomputedFile, job_dir),
        "spongeNetworkAnalysis": save_file(spongeNetworkAnalysis, job_dir),
        "spongeNetworkInteractions": save_file(
            spongeNetworkInteractions, job_dir
        ),
        "genieFile": save_file(genieFile, job_dir),
    }

    # 4) Build response payload
    payload = {
        "ok": True,
        "jobId": job_id,
        "email": str(email),
        "dataset": dataset.value,
        "method": method.value,
        "normalization": normalization,
        "filteringSpatial": filteringSpatial,
        "filteringSingleCell": filteringSingleCell,
        "tangram": tangram,
        "scores": scores_obj,
        "useDefaultLiana": useDefaultLiana,
        "options": options,
        "files": saved_files,
    }

    # 5) Persist a copy of the payload next to the uploaded files
    (job_dir / f"{job_id}_config.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )

    # 6) Return a clean JSON response the frontend can consume
    return payload


@app.post("/create_session/{name}")
async def create_session(name: str, response: Response):
    """
    Example: `curl -c cookies.txt -X POST http://127.0.0.1:3000/create_session/mopitas`
    """
    session = uuid4()
    data = SessionData(username=name)

    await backend.create(session, data)
    cookie.attach_to_response(response, session)

    return f"created session for {name}"


@app.get("/whoami", dependencies=[Depends(cookie)])
async def whoami(session_data: SessionData = Depends(verifier)):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:3000/whoami`
    """
    return session_data.username


@app.post("/delete_session")
async def del_session(response: Response, session_id: UUID = Depends(cookie)):
    """
    Example: `curl -b cookies.txt -X POST http://127.0.0.1:3000/delete_session`
    """
    await backend.delete(session_id)
    cookie.delete_from_response(response)
    return "deleted session"


@app.post("/read_adata")
async def read_adata(
    adata_path: AnnDataPath, session_id: UUID = Depends(cookie)
):
    """
    Example:
    ```
    curl -b cookies.txt \
    -H "Content-Type: application/json" \
    -X POST http://127.0.0.1:3000/read_adata \
    -d '{"path": "data/adata.h5ad"}'
    ```
    """
    session_data = await backend.read(session_id)

    adata = sc.read_h5ad(adata_path.path)

    reconstruct_obsm_cols = {
        "ligand_receptor_cosine_similarity": "ligand_receptor",
        "ligand_receptor_p_value": "ligand_receptor",
        "ligand_receptor_category": "ligand_receptor",
        "cell_comp_tf_activity_cosine_similarity": "cell_comp_tf_activity",
        "cell_comp_tf_activity_category": "cell_comp_tf_activity",
    }

    for obsm_key, col_names in reconstruct_obsm_cols.items():
        if obsm_key in adata.obsm:
            adata.obsm[obsm_key] = pd.DataFrame(
                adata.obsm[obsm_key],
                columns=adata.uns["liana_columns"][col_names],
                index=adata.obs_names,
            )

    session_data.adata = adata
    await backend.update(session_id, session_data)
    return {"status": "ok"}


@app.get("/obs/{column}", dependencies=[Depends(cookie)])
async def get_obs_column(
    column: str, session_data: SessionData = Depends(verifier)
):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:3000/obs/cell_type`
    """
    return session_data.adata.obs[column].to_dict()


@app.get("/var/{column}", dependencies=[Depends(cookie)])
async def get_var_column(
    column: str, session_data: SessionData = Depends(verifier)
):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:3000/var/n_cells`
    """
    return session_data.adata.var[column].to_dict()


@app.get("/obsm/{table}/{column}", dependencies=[Depends(cookie)])
async def get_obsm_column(
    table: str, column: str, session_data: SessionData = Depends(verifier)
):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:3000/obsm/ligand_receptor_cosine_similarity/LGALS9^PTPRC`
    """
    return session_data.adata.obsm[table][column].to_dict()


if __name__ == "__main__":
    app.state.data = None
    uvicorn.run(app, host="0.0.0.0", port=3000)
