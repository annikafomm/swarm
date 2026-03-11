# backend/main.py

from __future__ import annotations

# ---------------------------------
# Standard Library
# ---------------------------------
import json
import os
import shutil
import subprocess
import time
from dataset_management import DatasetRegistry
from contextlib import asynccontextmanager
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Optional
from uuid import UUID, uuid4

import numpy as np
import pandas as pd
import scanpy as sc
import uvicorn
from app import calculate_scores_helper

# merit
import re

# ---------------------------------
# Third-Party (FastAPI / Starlette)
# ---------------------------------
from fastapi import (
    Depends,
    FastAPI,
    File,
    Form,
    HTTPException,
    Request,
    Response,
    UploadFile,
)
from fastapi.middleware.cors import CORSMiddleware
from fastapi.exceptions import RequestValidationError
from fastapi.responses import FileResponse, JSONResponse
from fastapi_sessions.backends.implementations import InMemoryBackend
from fastapi_sessions.frontends.implementations import (
    CookieParameters,
    SessionCookie,
)
from fastapi_sessions.session_verifier import SessionVerifier
from pydantic import BaseModel as PydanticBaseModel
from pydantic import EmailStr
from starlette.responses import RedirectResponse
from datetime import datetime, timedelta
import asyncio

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

# Base folder for all uploads (created on startup).
# Use path relative to this file to avoid depending on current working dir.
#BASE_UPLOAD_DIR = Path(__file__).resolve().parent / "uploads"
BASE_UPLOAD_DIR = Path.cwd() / "../backend/uploads"
BASE_UPLOAD_DIR.mkdir(parents=True, exist_ok=True)

# Allow configuring CORS origins via environment variable (comma-separated).
# Example: ALLOWED_ORIGINS="https://myapp.com,https://staging.myapp.com"
ALLOWED_ORIGINS = [
    o.strip()
    for o in os.getenv(
        "ALLOWED_ORIGINS",
        "http://localhost:4200,http://127.0.0.1:4200,http://localhost:4201,http://127.0.0.1:4201"
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

class GeneSelectionMode(str, Enum):
    ctg = "ctg"
    hvg = "hvg"
    spapros = "spapros"
    svg = "svg"
    none = "None"

class Genome(str, Enum):
    Visium = "hg37"
    Xenium = "hg38"

class Method(str, Enum):
    Genie3 = "Genie3"
    Sponge = "Sponge"

class BaseModel(PydanticBaseModel):
    class Config:
        arbitrary_types_allowed = True


class SessionData(BaseModel):
    username: str
    adata_path: str = None
    genie_network_path: str | None = None
    sponge_network_path: str | None = None
    created_at: datetime = None

    def __init__(self, **data: Any):
        super().__init__(**data)
        if self.created_at is None:
            self.created_at = datetime.utcnow()



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


class SpongeNetworksPath(BaseModel):
    path: str


class GenieNetworkPath(BaseModel):
    path: str


# -----------------------------------------------------------------------------
# Application
# -----------------------------------------------------------------------------
# Lifespan event handler

@asynccontextmanager
async def lifespan(app: FastAPI):
    # Startup
    global dataset_registry
    dataset_registry = DatasetRegistry()

    # Register builtin datasets from backend/data (use absolute paths)
    base_path = Path(__file__).parent  # backend/ directory
    builtin_adata = base_path / "data" / "adata.h5ad"
    builtin_geojson = base_path.parent / "frontend" / "public" / "assets" / "hexagons.geojson"
    builtin_genie = base_path / "data" / "genie_network_filt.csv"
    builtin_sponge = base_path / "data" / "sponge_network_smaller.csv"

    print(f"Base path: {base_path}")
    print(f"Adata exists: {builtin_adata.exists()} at {builtin_adata}")

    if builtin_adata.exists():
        dataset_registry.register_builtin_dataset(
            dataset_id="builtin_main",
            alias="Default Dataset (Visium, BRCA)",
            adata_path=str(builtin_adata),
            geojson_path="/api/geojson/builtin_main",  # Use API URL for consistency
            genie_network_path=str(builtin_genie) if builtin_genie.exists() else None,
            sponge_network_path=str(builtin_sponge) if builtin_sponge.exists() else None,
            description="Pre-configured spatial transcriptomics dataset"
        )
        print(f"✓ Registered builtin dataset with paths:")
        print(f"  - adata: {builtin_adata}")
        print(f"  - geojson: /api/geojson/builtin_main")
        print(f"  - genie: {builtin_genie if builtin_genie.exists() else 'NOT FOUND'}")
        print(f"  - sponge: {builtin_sponge if builtin_sponge.exists() else 'NOT FOUND'}")
    else:
        print(f"✗ Builtin adata not found at {builtin_adata}")

    asyncio.create_task(cleanup_expired_sessions())
    yield
    # Shutdown
    pass

app = FastAPI(lifespan=lifespan)

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

def _load_adata_cached(file_path: str) -> sc.AnnData:
    """Load on-demand. Python's garbage collector will clean up when no longer referenced."""

    if file_path is None:
        # Fall back to builtin dataset if none is set
        base_path = Path(__file__).parent
        file_path = base_path / "data" / "adata.h5ad"
        if not file_path.exists():
            raise ValueError("No adata file has been loaded. Please call /read_adata first.")

    adata = sc.read_h5ad(str(file_path))
    reconstruct_obsm_cols = {
        "ligand_receptor_cosine_similarity": "ligand_receptor",
        "ligand_receptor_p_value": "ligand_receptor",
        "ligand_receptor_category": "ligand_receptor",
        "cell_comp_tf_activity_cosine_similarity": "cell_comp_tf_activity",
        "cell_comp_tf_activity_category": "cell_comp_tf_activity",
    }

    for obsm_key, col_names in reconstruct_obsm_cols.items():
        if obsm_key in adata.obsm and 'liana_columns' in adata.uns and col_names in adata.uns["liana_columns"]:
            adata.obsm[obsm_key] = pd.DataFrame(
                adata.obsm[obsm_key],
                columns=adata.uns["liana_columns"][col_names],
                index=adata.obs_names,
            )

    return adata

GRID_PREFIX = "grid_"

def strip_grid_prefix_dict(d: dict, prefix: str = GRID_PREFIX) -> dict:
    out = {}
    for k, v in d.items():
        if k.startswith(prefix):
            out[k[len(prefix):]] = v
        else:
            out[k] = v
    return out

def prefer_unprefixed_then_grid(adata, key: str, prefix: str = GRID_PREFIX) -> str | None:
    """
    Return the key that exists: prefer key, else prefix+key, else None
    """
    if key in adata.obs.columns:
        return key
    k2 = f"{prefix}{key}"
    if k2 in adata.obs.columns:
        return k2
    return None


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
    print(f"MAX_FILE_MB: {MAX_FILE_MB}, size_mb: {size_mb}")
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

    print(f"Saving uploaded file: {upload.filename}")
    _ensure_under_max_size(upload)
    print(f"ensured under max size: {upload.filename}")

    job_dir.mkdir(parents=True, exist_ok=True)
    safe_name = _sanitize_filename(upload.filename or "upload.bin")
    dest = job_dir / safe_name

    with dest.open("wb") as f:
        shutil.copyfileobj(upload.file, f)

    # Return as string for JSON serialization
    return str(dest.resolve())



def _step_and_borders_networks(edge_annotations: np.ndarray) -> dict[str, Any]:
    """
    Given an array of edge annotations, compute an appropriate step size and borders
    for binning the annotations into categories.

    Returns (step_size, min_border, max_border)
    """
    w_min = np.min(edge_annotations)
    w_max = np.max(edge_annotations)
    w_range = w_max - w_min

    if w_range == 0:
        return (1.0, w_min - 0.5, w_max + 0.5)

    # Determine step size based on range
    if w_range <= 0.1:
        step = 0.01
    elif w_range <= 1.0:
        step = 0.1
    elif w_range <= 10.0:
        step = 1.0
    else:
        step = 10.0

    num_possible_steps = w_range / step

    # Calculate borders
    min_border = np.floor(w_min / step) * step
    max_border = np.ceil(w_max / step) * step

    default_value = min_border + (np.ceil(num_possible_steps / 2) * step)

    return {"step": step, "min_border": min_border, "max_border": max_border, "default_value": default_value}



def get_subnetwork_data(file_path, gene_set, network_type):
    # gene_set is a set or list of gene names
    filtered_rows = []
    edge_annotations = []
    for chunk in pd.read_csv(file_path, chunksize=10000):
        if network_type == "genie":
            mask = chunk["regulatoryGene"].isin(gene_set) | chunk[
                "targetGene"
            ].isin(gene_set)
            annotation = "weight"
        elif network_type == "sponge":
            gene_mask = chunk["geneA"].isin(gene_set) | chunk["geneB"].isin(gene_set)

            mask = (
                gene_mask
                & chunk["p.adj"].notna()
                & chunk["p.adj"].lt(0.05)
                & chunk["mscor"].notna()
                & chunk["mscor"].abs().gt(0.1)
            )
            annotation = "p.adj"
        else:
            continue
        filtered_chunk = chunk[mask]
        edge_annotations.extend(filtered_chunk[annotation].values.tolist())
        if not filtered_chunk.empty:
            filtered_rows.append(filtered_chunk)
            # remove chunk from memory
            del filtered_chunk, chunk
    if filtered_rows:
        return pd.concat(filtered_rows, ignore_index=True), _step_and_borders_networks(np.array(edge_annotations))
    else:
        return pd.DataFrame()


# -----------------------------------------------------------------------------
# Convenience endpoints
# -----------------------------------------------------------------------------

async def cleanup_expired_sessions():
     while True:
        try:
            await asyncio.sleep(300)  # Check every 5 minutes
            current_time = datetime.now()
            expired_sessions = []

            # Iterate through all sessions
            for session_id, session_data in backend._data.items():
                age = current_time - session_data.created_at
                if age > timedelta(hours=1):  # Sessions expire after 1 hour
                    expired_sessions.append(session_id)

            # Delete expired sessions
            for session_id in expired_sessions:
                await backend.delete(session_id)
                print(f"Cleaned up expired session: {session_id}")
        except Exception as e:
            print(f"Error in cleanup task: {e}")



@app.post("/reset_session")
async def reset_session(session_id: UUID = Depends(cookie)):
    """Call this on page reload to clear old data"""
    session_data = await backend.read(session_id)
    # Reset all data except username
    session_data.adata_path = None
    session_data.genie_network_path = None
    session_data.sponge_network_path = None
    session_data.created_at = datetime.now()  # Reset timestamp
    await backend.update(session_id, session_data)
    return {"status": "session reset"}


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
@app.middleware("http")
async def log_requests(request: Request, call_next):
    print(f"[REQ] {request.method} {request.url} content-length={request.headers.get('content-length')}")
    return await call_next(request)

# Main endpoint for receiving form fields and files (multipart/form-data).
# It creates a dedicated job directory, saves all provided files there,
# stores a config JSON for reproducibility, and returns a summary payload.
@app.post("/api/upload", dependencies=[Depends(cookie)])
async def upload(
    # --- core ---
    email: Optional[str] = Form(None),
    dataset: Dataset = Form(...),
    # Spatial
    spatial_h5ad: UploadFile = File(...),
    spatial_normalization: bool = Form(False),
    spatial_filtering: bool = Form(False),
    # Tangram
    use_tangram: bool = Form(False),
    single_cell_h5ad: Optional[UploadFile] = File(None),
    singlecell_filtering: bool = Form(False),
    singlecell_normalization: bool = Form(False),
    gene_selection_mode: Optional[GeneSelectionMode] = Form(None),
    # Multiome
    use_multiome: bool = Form(False),
    multiome_rds: Optional[UploadFile] = File(None),

    # need_multiome_fragments: bool = Form(False),
    fragments_tsv_gz: Optional[UploadFile] = File(None),
    fragments_tsv_gz_tbi: Optional[UploadFile] = File(None),

    # Scores
    score_network: bool = Form(False),
    score_squidpy: bool = Form(False),
    score_liana_plus: bool = Form(False),
    score_chromVar: bool = Form(False),
    score_differential_motif_activity: bool = Form(False),
    score_motif_enrichment: bool = Form(False),
    score_FootprintingBias: bool = Form(False),

    genome: str = Form(None),

    # Network Algos
    alg_viper: bool = Form(False),
    alg_aucell: bool = Form(False),
    alg_gsva: bool = Form(False),
    alg_ssgsea: bool = Form(False),
    # SPONGE-Params
    net_m_score_threshold: Optional[float] = Form(None),
    net_p_adjust: Optional[str] = Form(None),
    net_ensembl_id_col: Optional[str] = Form(None),
    net_feature_col: Optional[str] = Form(None),
    net_rna_types: Optional[str] = Form(None),
    net_max_modules: Optional[int] = Form(None),
    # GENIE3-Params
    genie3_top_n_weights: Optional[int] = Form(None),
    genie3_n_regulatory_genes: Optional[int] = Form(None),
    genie3_n_regulons: Optional[int] = Form(None),
    # Network-Dateien
    genie3_network: Optional[UploadFile] = File(None),
    sponge_networkanalysis: Optional[UploadFile] = File(None),
    sponge_networkinteractions: Optional[UploadFile] = File(None),
    # Squidpy Flags + Params
    squidpy_moranI: bool = Form(False),
    squidpy_moranI_n_perms: Optional[int] = Form(None),
    squidpy_moranI_two_tailed: bool = Form(False),
    squidpy_moranI_corr_method: Optional[str] = Form(None),
    squidpy_gearyC: bool = Form(False),
    squidpy_gearyC_n_perms: Optional[int] = Form(None),
    squidpy_gearyC_two_tailed: bool = Form(False),
    squidpy_gearyC_corr_method: Optional[str] = Form(None),
    squidpy_centrality_score: bool = Form(False),
    squidpy_centrality_score_cluster_key: Optional[str] = Form(None),
    squidpy_co_occurrence: bool = Form(False),
    squidpy_co_occurrence_cluster_key: Optional[str] = Form(None),
    squidpy_co_occurrence_interval: Optional[int] = Form(None),
    squidpy_co_occurrence_n_splits: Optional[int] = Form(None),
    squidpy_neighborhood_enrichment: bool = Form(False),
    squidpy_neighborhood_enrichment_cluster_key: Optional[str] = Form(None),
    squidpy_neighborhood_enrichment_library_key: Optional[str] = Form(None),
    squidpy_neighborhood_enrichment_n_perms: Optional[int] = Form(None),
    # ChromVAR Flags + Params ---
    chromVar_moranI: bool = Form(False),
    chromVar_moranI_n_perms: Optional[int] = Form(None),
    chromVar_moranI_two_tailed: Optional[str] = Form(None),   # "oneTailed" / "twoTailed"
    chromVar_moranI_corr_method: Optional[str] = Form(None),

    chromVar_gearyC: bool = Form(False),
    chromVar_gearyC_n_perms: Optional[int] = Form(None),
    chromVar_gearyC_two_tailed: Optional[str] = Form(None),
    chromVar_gearyC_corr_method: Optional[str] = Form(None),

    chromVar_differential_motif_activity: bool = Form(False),

    # LIANA
    liana_composition_column: Optional[str] = Form(None),
    liana_genie3_network: Optional[UploadFile] = File(None),
    liana_pathway_network: Optional[UploadFile] = File(None),
    session_data: "SessionData" = Depends(verifier),
):
    print("in method")
    raw_username = session_data.username
    #raw_username = "merit"
    user_safe = _sanitize_filename(raw_username) or "anon"
    job_id = f"job_{int(time.time() * 1000)}"
    job_dir = BASE_UPLOAD_DIR / f"{job_id}_{user_safe}"
    job_dir.mkdir(parents=True, exist_ok=True)

    with open(job_dir / f"{job_id}_multiome.txt", "w") as f:
        f.write(f"use_multiome: {use_multiome}\n")

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

    # 3) Save all files into that directory
    # Dateien speichern (alles in ./uploads)
    saved_files = {
        "spatial_h5ad": save_file(spatial_h5ad, job_dir),
        "single_cell_h5ad": save_file(single_cell_h5ad, job_dir),
        "multiome_rds": save_file(multiome_rds, job_dir),
        "fragments_tsv_gz": save_file(fragments_tsv_gz, job_dir),
        "fragments_tsv_gz_tbi": save_file(fragments_tsv_gz_tbi, job_dir),
        "genie3_network": save_file(genie3_network, job_dir),
        "sponge_networkanalysis": save_file(sponge_networkanalysis, job_dir),
        "sponge_networkinteractions": save_file(
            sponge_networkinteractions, job_dir
        ),
        "liana_genie3_network": save_file(liana_genie3_network, job_dir),
        "liana_pathway_network": save_file(liana_pathway_network, job_dir),
    }


    # TODO: if .rds: convert to h5ad file, and read again
    if use_multiome:
        print("in rds2h5ad block")
        rds_path = saved_files["multiome_rds"]
        assay_name = "RNA"
        h5ad_path = re.sub(r"\.rds$", f".h5ad", rds_path)
        print(f"Running file conversion in multiome env ...")
        log_path = Path(h5ad_path).with_suffix(".log")
        with log_path.open("w") as log_file:
            result = subprocess.run(
                [   "Rscript",
                    "../backend/rds_to_h5ad.R",
                    "--rds_path", rds_path,
                    "--assay", assay_name,
                    "--h5ad_path", h5ad_path,
                ],
                stdout=log_file,
                stderr=log_file,
                text=True,
                check=False,
            )

        use_tangram = True
        single_cell_h5ad = h5ad_path
        print(single_cell_h5ad)
        saved_files["single_cell_h5ad"] = single_cell_h5ad


    # TODO: change saved files dict. spatial_h5ad should become the h5ad path
    # or single_cell_h5ad?

    # 4) Build response payload
    payload = {
        "email": str(email),
        "dataset": dataset.value,
        "spatial": {
            "normalization": spatial_normalization,
            "filtering": spatial_filtering,
        },
        "files": saved_files,
        "tangram": {
            "use": use_tangram,
            "filtering": singlecell_filtering if use_tangram else None,
            "normalization": singlecell_normalization if use_tangram else None,
            "gene_selection_mode": (
                gene_selection_mode.value if (use_tangram and gene_selection_mode) else None
            ),
        },
        "multiome": {
            "use": use_multiome,
        },
        "scores": {
            "network": score_network,
            "squidpy": score_squidpy,
            "liana_plus": score_liana_plus,
            "chromVar": score_chromVar,
            "differential_motif_activity": score_differential_motif_activity,
            "motif_enrichment": score_motif_enrichment,
            "footprinting": score_FootprintingBias,
        },
        "genome": genome,

        "network": {
            "algorithms": {
                "viper": alg_viper,
                "aucell": alg_aucell,
                "gsva": alg_gsva,
                "ssgsea": alg_ssgsea,
            },
            "sponge_params": {
                "m_score_threshold": net_m_score_threshold,
                "p_adjust": net_p_adjust,
                "ensembl_id_col": net_ensembl_id_col,
                "feature_col": net_feature_col,
                "rna_types": net_rna_types,
                "max_modules": net_max_modules,
            },
            "genie3_params": {
                "top_n_weights": genie3_top_n_weights,
                "n_regulatory_genes": genie3_n_regulatory_genes,
                "n_regulons": genie3_n_regulons,
            },
        },
        "squidpy": {
            "moranI": squidpy_moranI,
            "moranI_params": {
                "n_perms": squidpy_moranI_n_perms,
                "two_tailed": squidpy_moranI_two_tailed,
                "corr_method": squidpy_moranI_corr_method,
            },
            "gearyC": squidpy_gearyC,
            "gearyC_params": {
                "n_perms": squidpy_gearyC_n_perms,
                "two_tailed": squidpy_gearyC_two_tailed,
                "corr_method": squidpy_gearyC_corr_method,
            },
            "centrality_score": squidpy_centrality_score,
            "centrality_score_params": {
                "cluster_key": squidpy_centrality_score_cluster_key
            },
            "co_occurrence": squidpy_co_occurrence,
            "co_occurrence_params": {
                "cluster_key": squidpy_co_occurrence_cluster_key,
                "interval": squidpy_co_occurrence_interval,
                "n_splits": squidpy_co_occurrence_n_splits,
            },
            "neighborhood_enrichment": squidpy_neighborhood_enrichment,
            "neighborhood_enrichment_params": {
                "cluster_key": squidpy_neighborhood_enrichment_cluster_key,
                "library_key": squidpy_neighborhood_enrichment_library_key,
                "n_perms": squidpy_neighborhood_enrichment_n_perms,
            },
        },
        "chromVar": {
            "moranI": chromVar_moranI,
            "moranI_params": {
                "n_perms": chromVar_moranI_n_perms,
                # bool in downstream code if you like:
                "two_tailed": (chromVar_moranI_two_tailed == "twoTailed"),
                "tails": chromVar_moranI_two_tailed,  # keep raw if useful
                "corr_method": chromVar_moranI_corr_method,
            },
            "gearyC": chromVar_gearyC,
            "gearyC_params": {
                "n_perms": chromVar_gearyC_n_perms,
                "two_tailed": (chromVar_gearyC_two_tailed == "twoTailed"),
                "tails": chromVar_gearyC_two_tailed,
                "corr_method": chromVar_gearyC_corr_method,
            },
            "differential_motif_activity": chromVar_differential_motif_activity,
        },
        "liana": {
            "composition_column": liana_composition_column,
        },
    }

    # TODO call visium_to_geojson

    out_dir = await calculate_scores_helper(job_dir, payload)

    # Xenium: scores are computed on a grid/spot representation, then broadcast back to cells.
    # For the web app we prefer the *cell-level* AnnData (xenium_cells_with_grid_scores.h5ad).
    # If not available, we fall back to spot-level (st_scores.h5ad) or the original upload.


    print("\n=== DEBUG upload() ===")
    print("dataset:", dataset.value)
    print("job_dir:", job_dir)
    print("out_dir:", out_dir)

    if out_dir is not None:
        print("files in out_dir:", sorted(os.listdir(out_dir)))
    else:
        print("out_dir is None (pipeline returned None)")

    adata_path = None
    selected_reason = None
    tangram_adata_path = None
    out_files = {}

    if out_dir is not None:
        for filename in os.listdir(out_dir):
            if filename.endswith("xenium_cells_with_grid_scores.h5ad") and dataset == Dataset.Xenium:
                adata_path = os.path.join(out_dir, filename)
                selected_reason = "Found xenium cell-level output"
                break

        if adata_path is None:
            for filename in os.listdir(out_dir):
                if filename.endswith("st_scores.h5ad"):
                    adata_path = os.path.join(out_dir, filename)
                    selected_reason = "Fallback to grid-level st_scores (cell-level missing)"
                    break

        if filename.endswith("tg_scores.h5ad"):
            tangram_adata_path = os.path.join(out_dir, filename)

        if filename.endswith("xenium_tg_scores.h5ad"):
            tangram_adata_path = os.path.join(out_dir, filename)

        print("selected adata_path:", adata_path)
        print("selection reason:", selected_reason)




    if adata_path is None:
        adata_path = saved_files.get("spatial_h5ad")
        selected_reason = "Fallback to uploaded spatial_h5ad (no pipeline output found)"
        print("WARNING:", selected_reason)

    if out_dir is None:
        out_dir = job_dir

    if adata_path is not None:
        dataset_id = f"{job_id}_{user_safe}"

        email_prefix = str(email.split("@")[0]) if email else "unknown"
        timestamp = datetime.now().strftime("%Y%m%d_%H%M%S")

        dataset_alias = f"{email_prefix}_{dataset_id}_{timestamp}"
        out_files["adataPath"] = adata_path

        data_type = dataset.value.lower()
        footprint_list = [
            str(Path(out_dir).relative_to(BASE_UPLOAD_DIR) / f)
            for f in os.listdir(out_dir)
            if f.startswith("footprint") and f.endswith(".pdf")
        ]
        dataset_registry.register_uploaded_dataset(
            dataset_id=dataset_id,
            alias=dataset_alias,
            adata_path=adata_path,
            tangram_adata_path=tangram_adata_path if use_tangram else None,
            genie_network_path=out_files.get("genie_network_path"),
            sponge_network_path=out_files.get("sponge_network_path"),
            geojson_path=f"/api/geojson/{dataset_id}",  # Use API URL instead of file path
            user=user_safe,
            dataset_type=dataset,
            use_tangram=use_tangram,
            created_at=datetime.now().isoformat(),
            footprint_list=footprint_list,
            # TODO Add paths for multiome
        )

        geojson_out = os.path.join(job_dir, "hexagons.geojson")
        print("geojson_out:", geojson_out)
        print("geojson input adata:", adata_path)
        print("geojson data_type:", dataset.value.lower())

        # Convert selected AnnData to GeoJSON for map rendering.
        subprocess.run(
            [
                "python",
                "../backend/visium_to_geojson.py",
                "--adata",
                adata_path,
                "--outpath",
                geojson_out,
                "--data_type",
                data_type,
            ]
        )
        out_files["geojsonPath"] = os.path.join(out_dir, "hexagons.geojson")
        print("geojson exists after run?:", os.path.isfile(geojson_out))

    if os.path.isfile(os.path.join(out_dir, "genie_network_filtered_st.csv")):
        out_files["genieFiltPath"] = os.path.join(out_dir, "genie_network_filtered_st.csv")

    if os.path.isfile(os.path.join(out_dir, "sponge_network_filtered_st.csv")):
        out_files["spongeFiltPath"] = os.path.join(out_dir, "sponge_network_filtered_st.csv")

    payload["output_files"] = out_files


    # 5) Persist a copy of the payload next to the uploaded files
    (job_dir / f"{job_id}_config.json").write_text(
        json.dumps(payload, indent=2), encoding="utf-8"
    )

    # 6) Return a clean JSON response the frontend can consume
    return payload



@app.get("/api/datasets", dependencies=[Depends(cookie)])
async def get_datasets(session_data: SessionData = Depends(verifier)):
    """
    Return all available datasets for the current user + builtin datasets + shared datasets.
    Filters out datasets with missing files.
    """
    try:
        all_datasets = dataset_registry.get_all_datasets()
        print(f"All datasets: {all_datasets}")
        user_datasets = dataset_registry.get_user_datasets(session_data.username)
        print(f"User datasets for {session_data.username}: {user_datasets}")

        # Verify files exist before returning
        valid_user_datasets = {}
        for dataset_id, info in user_datasets.items():
            if Path(info["adata_path"]).exists():
                valid_user_datasets[dataset_id] = info

        # Also include shared datasets (marked with user="__shared__")
        shared_datasets = {}
        if "uploaded" in all_datasets:
            for dataset_id, info in all_datasets["uploaded"].items():
                if info.get("user") == "__shared__" and Path(info["adata_path"]).exists():
                    shared_datasets[dataset_id] = info

        datasets_json = {
            "builtin": all_datasets.get("builtin", {}),
            "uploaded": {**valid_user_datasets, **shared_datasets}
        }

        print(f"Returning datasets: {datasets_json}")
        return datasets_json
    except Exception as e:
        print(f"✗ Error in get_datasets: {e}")
        import traceback
        traceback.print_exc()
        raise


@app.get("/api/geojson/{dataset_id}", dependencies=[Depends(cookie)])
async def get_geojson(dataset_id: str):
    """Serve GeoJSON files for datasets"""
    try:
        # Construct path from dataset_id - uploaded datasets have format: {job_id}_{user_safe}
        # Builtin datasets have format: builtin_{name}

        if dataset_id.startswith("builtin_"):
            # Builtin dataset - located in frontend/public/assets/
            geojson_path = Path(__file__).parent.parent / "frontend" / "public" / "assets" / "hexagons.geojson"
            print(f"[DEBUG] Looking for builtin geojson at: {geojson_path}")
        else:
            # Uploaded dataset - extract job_id from dataset_id format: job_TIMESTAMP_USER
            # The directory is created as: job_TIMESTAMP_USER
            geojson_path = BASE_UPLOAD_DIR / dataset_id / "hexagons.geojson"
            print(f"[DEBUG] Looking for uploaded geojson at: {geojson_path}")

        print(f"[DEBUG] File exists: {geojson_path.exists()}")

        if not geojson_path.exists():
            print(f"[ERROR] GeoJSON file not found at {geojson_path}")
            raise HTTPException(status_code=404, detail=f"GeoJSON file not found: {geojson_path}")

        print(f"[DEBUG] Serving geojson from: {geojson_path}")
        return FileResponse(
            geojson_path,
            media_type="application/geo+json",
            filename=f"{dataset_id}.geojson"
        )
    except HTTPException:
        raise
    except Exception as e:
        print(f"[ERROR] Exception in get_geojson: {e}")
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))
async def get_hexagon(user: str, subdir: str, filename: str):
    """
    file_path = Path("./uploads") / filename
    if file_path.exists() and file_path.is_file():
        return FileResponse(str(file_path))
    """
    file_path = BASE_UPLOAD_DIR / Path(user) / Path(subdir) / filename
    if file_path.exists() and file_path.is_file():
        return FileResponse(str(file_path))
    raise HTTPException(status_code=404, detail="File not found")


# ============================================================================
# DEBUG / TESTING ENDPOINTS (Backend-only, not exposed to frontend)
# ============================================================================

@app.post("/debug/rescan_uploads")
async def debug_rescan_uploads():
    """
    TESTING ONLY: Re-scan the uploads folder and re-register existing datasets.
    Useful when you've run the pipeline and want to reload without restarting the server.

    Example: `curl -X POST http://127.0.0.1:3000/debug/rescan_uploads`

    Note: This endpoint is intentionally backend-only and not exposed in the frontend.
    """
    try:
        # Rescan the uploads folder
        results = dataset_registry.rescan_uploads_folder(BASE_UPLOAD_DIR)

        # Reload datasets into memory
        # Re-load to get the latest from disk
        dataset_registry.datasets = dataset_registry._load_registry()

        return {
            "status": "success",
            "message": "Rescanned uploads folder and reregistered datasets",
            "results": results,
            "total_datasets": len(dataset_registry.get_all_datasets().get("uploaded", {}))
        }
    except Exception as e:
        print(f"Error in debug rescan: {e}")
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))


@app.post("/debug/rescan_uploads_shared")
async def debug_rescan_uploads_shared():
    """
    TESTING ONLY: Rescan uploads folder and register ALL datasets as shared.
    Datasets marked as shared will appear in the dropdown for all users.

    This is useful for testing when you want datasets to be accessible without
    worrying about session/user filtering.

    Example: `curl -X POST http://127.0.0.1:3000/debug/rescan_uploads_shared`
    """
    try:
        # Rescan the uploads folder
        results = dataset_registry.rescan_uploads_folder(BASE_UPLOAD_DIR)

        # Now convert all newly registered datasets to "shared"
        all_datasets = dataset_registry.get_all_datasets()
        if "uploaded" in all_datasets:
            for dataset_id, info in all_datasets["uploaded"].items():
                # Change user to __shared__ so it appears for everyone
                info["user"] = "__shared__"

        # Save the updated registry
        dataset_registry._save_registry()

        # Reload datasets into memory without clearing shared datasets
        dataset_registry.datasets = dataset_registry._load_registry(clear_uploads=False)

        return {
            "status": "success",
            "message": "Rescanned uploads and marked all datasets as SHARED (visible to all users)",
            "results": results,
            "total_shared_datasets": len(dataset_registry.get_all_datasets().get("uploaded", {}))
        }
    except Exception as e:
        print(f"Error in debug rescan shared: {e}")
        import traceback
        traceback.print_exc()
        raise HTTPException(status_code=500, detail=str(e))


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
    try:
        await backend.delete(session_id)
    except KeyError:
        pass

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

    session_data.adata_path = adata_path.path

    print(f"Setting adata path to: {adata_path.path}")

    await backend.update(session_id, session_data)
    return {"status": "ok"}


@app.get("/download_adata", dependencies=[Depends(cookie)])
async def download_adata(session_data: SessionData = Depends(verifier)):
    file_path = session_data.adata_path
    if not os.path.exists(file_path):
        return {"error": "File not found"}
    return FileResponse(
        file_path,
        media_type="application/octet-stream",
        filename="map.h5ad",
    )


@app.post("/read_network_genie")
async def read_network_genie(
    genie_network_path: GenieNetworkPath, session_id: UUID = Depends(cookie)
):
    """
    Example:
    ```
    curl -b cookies.txt \
    -H "Content-Type: application/json" \
    -X POST http://127.0.0.1:3000/read_networks \
    -d '{"path": "data/networks.json"}'
    ```
    """

    session_data = await backend.read(session_id)

    session_data.genie_network_path = genie_network_path.path
    await backend.update(session_id, session_data)

    return {"status": "ok"}


@app.post("/read_network_sponge")
async def read_network_sponge(sponge_network_path: SpongeNetworksPath, session_id: UUID = Depends(cookie)
):
    """
    Example:
    ```
    curl -b cookies.txt \
    -H "Content-Type: application/json" \
    -X POST http://127.0.0.1:3000/read_sponge_network \
    -d '{"path": "data/sponge_network.json"}'
    ```
     """
    session_data = await backend.read(session_id)

    session_data.sponge_network_path = sponge_network_path.path

    await backend.update(session_id, session_data)

    return {"status": "ok"}


@app.get("/geneset_connections_genie", dependencies=[Depends(cookie)])
async def get_geneset_connections(
    gene_set_name: str, session_data: SessionData = Depends(verifier)
):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:3000/geneset_connections`
    """

    if not session_data.genie_network_path:
        return {"connections": [], "slider_data": {}}


    # Get the geneset from the name
    adata = _load_adata_cached(session_data.adata_path)

    # Check if genie_genesets exists in adata.uns
    if "genie_genesets" not in adata.uns:
        return {"connections": [], "slider_data": {}}

    gene_set = adata.uns["genie_genesets"].get(
        gene_set_name, None
    )
    gene_set = list(gene_set) + [gene_set_name]

    # Get connections from genie_network
    connections, slider_data = get_subnetwork_data(
        session_data.genie_network_path, gene_set, "genie"
    )

    # Write to dict
    connections = connections.to_dict(orient="records")

    return {"connections": connections, "slider_data": slider_data}


@app.get("/geneset_connections_sponge", dependencies=[Depends(cookie)])
async def get_geneset_connections(
    gene_set_name: str, session_data: SessionData = Depends(verifier)
):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:3000/geneset_connections`
    """
    # Get the geneset from the name
    if not session_data.sponge_network_path:
        return {"connections": [], "slider_data": {}}

    adata = _load_adata_cached(session_data.adata_path)
    gene_set = adata.uns["sponge_genesets"].get(
        gene_set_name, None
    )
    gene_set = list(gene_set) + [gene_set_name]

    # Get connections from sponge_network
    connections, slider_data = get_subnetwork_data(
        session_data.sponge_network_path, gene_set, "sponge"
    )
    connections = connections.to_dict(orient="records")

    return {"connections": connections, "slider_data": slider_data}


@app.get("/obs/{column}", dependencies=[Depends(cookie)])
async def get_obs_column(
    column: str, session_data: SessionData = Depends(verifier)
):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:3000/obs/cell_type`
    """
    adata = _load_adata_cached(session_data.adata_path)
    return adata.obs[column].to_dict()



@app.get("/var/{column}", dependencies=[Depends(cookie)])
async def get_var_column(
    column: str, session_data: SessionData = Depends(verifier)
):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:3000/var/n_cells`
    """
    adata = _load_adata_cached(session_data.adata_path)
    return adata.var[column].to_dict()


# @app.get("/obsm/{table}/{column}", dependencies=[Depends(cookie)])
# async def get_obsm_column(
#     table: str, column: str, session_data: SessionData = Depends(verifier)
# ):
#     """
#     Example: `curl -b cookies.txt http://127.0.0.1:3000/obsm/ligand_receptor_cosine_similarity/LGALS9^PTPRC`
#     """
#     obsm_data = session_data.adata.obsm[table]

#     # Special handling for chromvar_spot_scores: convert to DataFrame with motif names
#     if table == "chromvar_spot_scores" and "chromvar_motifs" in session_data.adata.uns:
#         motif_names = list(session_data.adata.uns["chromvar_motifs"])
#         if not isinstance(obsm_data, pd.DataFrame):
#             obsm_data = pd.DataFrame(
#                 obsm_data,
#                 index=session_data.adata.obs_names,
#                 columns=motif_names
#             )

#     return obsm_data[column].to_dict()

@app.get("/obsm/{table}/{column}", dependencies=[Depends(cookie)])
async def get_obsm_column(table: str, column: str, session_data: SessionData = Depends(verifier)):
    print(f"[DEBUG] Endpoint called: table={table}, column={column}")
    print(f"[DEBUG] Session username: {session_data.username}")
    adata = _load_adata_cached(session_data.adata_path)

    # --- ChromVAR special case: column is motif name OR comma-separated motif names ---
    if table == "chromvar_spot_scores":
        try:
            if adata is None:
                raise HTTPException(status_code=500, detail="No adata loaded in session")

            if "chromvar_motifs" not in adata.uns:
                print(f"[DEBUG] Available uns keys: {list(adata.uns.keys())}")
                raise HTTPException(status_code=500, detail="adata.uns['chromvar_motifs'] missing")

            chromvar_motifs = np.asarray(adata.uns["chromvar_motifs"])
            print(f"[DEBUG] chromvar_motifs shape: {chromvar_motifs.shape}, dtype: {chromvar_motifs.dtype}")
            print(f"[DEBUG] First 5 motifs: {chromvar_motifs[:5]}")

            motif_list = column.split(",")  # supports "MA1,MA2,MA3"
            print(f"[DEBUG] Requested motifs: {motif_list}")

            idx = np.where(np.isin(chromvar_motifs, motif_list))[0]
            print(f"[DEBUG] Found indices: {idx}, size: {idx.size}")

            if idx.size == 0:
                raise HTTPException(status_code=404, detail=f"No motifs found from: {motif_list}")

            if "chromvar_spot_scores" not in adata.obsm:
                print(f"[DEBUG] Available obsm keys: {list(adata.obsm.keys())}")
                raise HTTPException(status_code=500, detail="chromvar_spot_scores not in adata.obsm")

            X = adata.obsm["chromvar_spot_scores"]
            print(f"[DEBUG] chromvar_spot_scores type: {type(X)}, shape: {X.shape if hasattr(X, 'shape') else 'no shape'}")

            # sum scores across selected motifs
            scores = np.asarray(X[:, idx].sum(axis=1)).ravel()
            print(f"[DEBUG] scores shape: {scores.shape}, dtype: {scores.dtype}")

            return {bc: float(s) for bc, s in zip(adata.obs_names, scores)}
        except HTTPException:
            raise
        except Exception as e:
            print(f"[ERROR] chromvar_spot_scores endpoint failed: {type(e).__name__}: {e}")
            import traceback
            traceback.print_exc()
            raise HTTPException(status_code=500, detail=f"ChromVAR processing failed: {str(e)}")

    # --- default behavior unchanged ---
    obsm_data = adata.obsm[table]
    if not isinstance(obsm_data, pd.DataFrame):
        obsm_data = pd.DataFrame(obsm_data, index=adata.obs_names)

    return obsm_data[column].to_dict()

@app.get("/obsm/regulatory_scores/cell/{barcode}", dependencies=[Depends(cookie)])
async def get_obsm_row(
    barcode: str, session_data: SessionData = Depends(verifier)
):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:3000/obsm/AAACCTCATGAAGTTG-1`
    """
    regulatory_scores = [
    'aucell_scores',
    'spongeeeffects_ssGSEA_scores',
    'spongeeffects_GSVA_scores',
    'viper_scores',
    ]
    adata = _load_adata_cached(session_data.adata_path)
    available_scores = adata.obsm.keys()
    available_scores = [score for score in available_scores if score.endswith("_genie3") or score.endswith("_sponge")]

    row_data = {}
    for score in available_scores:
        obsm_data = adata.obsm[score]
        if not isinstance(obsm_data, pd.DataFrame):
            obsm_data = pd.DataFrame(
                obsm_data, index=adata.obs_names
            )
        row_data[score] = obsm_data.loc[barcode].to_dict()
    return row_data



@app.get("/api/download/{file_path:path}", dependencies=[Depends(cookie)])
async def download_file(file_path: str):
    full_path = (BASE_UPLOAD_DIR / file_path).resolve()

    # TODO: remove this. just so that hardcodded footprint plot can be used
    # Prevent path traversal outside uploads dir
    # if BASE_UPLOAD_DIR not in full_path.parents:
    #     raise HTTPException(status_code=403, detail="Forbidden")

    if not full_path.exists() or not full_path.is_file():
        raise HTTPException(status_code=404, detail=f"File not found: {file_path}")

    return FileResponse(str(full_path))




# @app.get("/api/footprints/{user}/{subdir}")
# async def list_footprints(user: str, subdir: str):
#     job_dir = BASE_UPLOAD_DIR / Path(user) / Path(subdir)

#     if not job_dir.exists() or not job_dir.is_dir():
#         raise HTTPException(status_code=404, detail="Directory not found")

#     files = sorted(
#         f.name for f in job_dir.glob("footprint_*.pdf") if f.is_file()
#     )
#     return files


# @app.get("/api/footprints/{user}/{subdir}/{filename}")
# async def get_footprint(user: str, subdir: str, filename: str):

#     safe_name = Path(filename).name
#     if safe_name != filename:
#         raise HTTPException(status_code=400, detail="Invalid filename")

#     if not safe_name.startswith("footprint_") or not safe_name.endswith(".pdf"):
#         raise HTTPException(status_code=404, detail="File not found")

#     file_path = BASE_UPLOAD_DIR / Path(user) / Path(subdir) / safe_name

#     if file_path.exists() and file_path.is_file():
#         return FileResponse(str(file_path), media_type="application/pdf")

#     raise HTTPException(status_code=404, detail="File not found")


@app.get("/X/{gene}", dependencies=[Depends(cookie)])
async def get_X_by_gene(
    gene: str, session_data: SessionData = Depends(verifier)
):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:3000/X/ENSG00000241860`
    """
    adata = _load_adata_cached(session_data.adata_path)
    expressions = adata[:, gene].X.toarray().flatten().tolist()
    barcodes = adata.obs.index
    return {
        barcode: expression
        for barcode, expression in zip(barcodes, expressions)
    }


# ---------------------------------------------------------------------------
# On-demand footprint computation endpoints
# ---------------------------------------------------------------------------

@app.get("/api/motifs", dependencies=[Depends(cookie)])
async def get_available_motifs(
    session_data: SessionData = Depends(verifier),
    dataset_id: Optional[str] = None,
):
    """
    Return the list of chromVAR motif IDs for the current dataset.

    Resolution order:
    1. Active session (session_data.adata_path set via /read_adata) — normal flow.
    2. Fallback: dataset_id query param → look up adata_path from the registry.
       This handles the rescanned-dataset case where /read_adata was never called.
    """
    adata_path = session_data.adata_path

    if not adata_path and dataset_id:
        ds = dataset_registry.get_dataset_by_id(dataset_id)
        if not ds:
            raise HTTPException(status_code=404, detail=f"Dataset '{dataset_id}' not found in registry")
        adata_path = ds.get("adata_path")

    if not adata_path:
        raise HTTPException(
            status_code=400,
            detail="No dataset loaded in session and no dataset_id provided. "
                   "Either call /read_adata or pass ?dataset_id=<id>."
        )

    try:
        adata = _load_adata_cached(adata_path)
        motifs = list(adata.uns.get("chromvar_motifs", []))
        print(f"[get_available_motifs] Found motifs[0:5]: {motifs[:5]}, total={len(motifs)}, dataset_id={dataset_id}")
        return {"motifs": motifs}
    except Exception as e:
        raise HTTPException(status_code=500, detail=f"Failed to load motifs: {e}")


@app.post("/api/compute_footprint", dependencies=[Depends(cookie)])
async def compute_footprint(
    motif: str = Form(...),
    cluster_by: str = Form("cell_type"),
    dataset_id: Optional[str] = Form(None),
    session_data: SessionData = Depends(verifier),
):
    """
    Trigger on-demand footprint computation for a specific motif.
    Runs compute_additional_footprints.R and returns the new footprint URL.

    Resolution order for adata_path:
    1. Active session (session_data.adata_path set via /read_adata).
    2. Fallback: dataset_id form field → registry lookup (rescanned-dataset case).
    """
    adata_path = session_data.adata_path

    if not adata_path and dataset_id:
        ds = dataset_registry.get_dataset_by_id(dataset_id)
        if not ds:
            raise HTTPException(status_code=404, detail=f"Dataset '{dataset_id}' not found in registry")
        adata_path = ds.get("adata_path")

    if not adata_path:
        raise HTTPException(
            status_code=400,
            detail="No dataset loaded in session and no dataset_id provided. "
                   "Either call /read_adata or pass dataset_id in the form."
        )

    # The job output dir is the directory containing the adata file
    out_dir = Path(adata_path).resolve().parent

    # Validate required files exist before launching R
    spot_rds = out_dir / "spot_obj_footprints.rds"
    dissociated_rds = out_dir / "dissociated_obj_footprints.rds"
    if not spot_rds.exists() or not dissociated_rds.exists():
        raise HTTPException(
            status_code=422,
            detail=(
                "Pre-computed footprint RDS objects not found in the job directory. "
                "The initial pipeline must be run with --footprinting first."
            ),
        )


    import subprocess
    cmd = [
        "Rscript", "../backend/calc_multiome_scores/compute_additional_footprints.R",
        "--outdir", str(out_dir),
        "--motif", motif,
        "--cluster_by", cluster_by,
    ]
    print(f"[compute_footprint] Running: {' '.join(cmd)}")
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,  # merge stderr into stdout so both stream together
        text=True,
        bufsize=1,  # line-buffered
    )
    output_lines: list[str] = []
    for line in proc.stdout:  # type: ignore[union-attr]
        print(f"[R] {line}", end="", flush=True)
        output_lines.append(line)
    proc.wait()
    if proc.returncode != 0:
        full_output = "".join(output_lines)
        raise HTTPException(
            status_code=500,
            detail=f"R script failed (exit {proc.returncode}): {full_output[-2000:]}"
        )

    # Verify the PDF was actually produced
    pdf_filename = f"footprint_{motif}.pdf"
    pdf_path = out_dir / pdf_filename
    if not pdf_path.exists():
        raise HTTPException(status_code=500, detail=f"Expected output {pdf_filename} was not created")

    # TODO: decide whether to add here
    # Update the dataset registry entry that owns this adata file
    relative_pdf = str(pdf_path.relative_to(BASE_UPLOAD_DIR.resolve()))
    all_datasets = dataset_registry.get_all_datasets()
    for category in ("uploaded", "builtin"):
        for ds_id, info in all_datasets.get(category, {}).items():
            if info.get("adata_path") == adata_path:
                existing = info.get("footprint_list") or []
                if relative_pdf not in existing:
                    dataset_registry.update_dataset_paths(ds_id, footprint_list=existing + [relative_pdf])
                break

    return {
        "footprint_url": f"/api/download/{relative_pdf}",
        "relative_path": relative_pdf,
        "motif": motif,
        "cluster_by": cluster_by,
    }


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=3000)
    # for merit
    #uvicorn.run(app, host="0.0.0.0", port=3005)
