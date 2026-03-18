# backend/main.py

import json
import os
import shutil
import subprocess
import time
from dataset_management import DatasetRegistry
from dataset_structure import Params, DatasetFactory
from contextlib import asynccontextmanager
from enum import Enum
from pathlib import Path
from typing import Any, Dict, List, Optional
from uuid import UUID, uuid4

import numpy as np
import pandas as pd
import scanpy as sc
import uvicorn
from app import calculate_scores_helper
from models import (
    UploadRequest,
    UploadResponse,
    OutputFiles,
    FilesInput,
    DatasetType,
    GeneSelectionMode,
)

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
from scipy import sparse

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
MAX_FILE_MB: Optional[int] = 5000

# TTL cache for per-(dataset,gene) expression min/max stats.
# Key: (dataset_id, normalized_gene)
# Value: {"min": float, "max": float, "expires_at": datetime}
GENE_STATS_CACHE: Dict[tuple[str, str], Dict[str, Any]] = {}
GENE_STATS_CACHE_TTL_SECONDS = 600

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


# Simple models for path parameters used in specific endpoints
class AnnDataPath(BaseModel):
    path: str


class SpongeNetworksPath(BaseModel):
    path: str


class GenieNetworkPath(BaseModel):
    path: str


# =============================================================================
# Application
# =============================================================================
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


def _normalize_gene_key(gene: str) -> str:
    return gene.strip().lower()


def _get_visible_dataset_paths(session_data: SessionData) -> Dict[str, str]:
    """
    Return dataset_id -> adata_path for datasets visible to this user.
    Visibility:
      - builtin datasets (always)
      - uploaded datasets owned by user or marked shared (__shared__)
    """
    visible: Dict[str, str] = {}
    all_datasets = dataset_registry.get_all_datasets(as_dict=True)

    for dataset_id, dataset_dict in all_datasets.get("builtin", {}).items():
        adata_path = dataset_dict.get("adata_path")
        if adata_path and Path(adata_path).exists():
            visible[dataset_id] = adata_path

    username = session_data.username
    for dataset_id, dataset_dict in all_datasets.get("uploaded", {}).items():
        owner = dataset_dict.get("user")
        if owner not in {username, "__shared__"}:
            continue
        adata_path = dataset_dict.get("adata_path")
        if adata_path and Path(adata_path).exists():
            visible[dataset_id] = adata_path

    return visible


def _resolve_adata_path(session_data: SessionData, dataset_id: Optional[str] = None) -> str:
    """
    Resolve adata path with optional dataset override.

    Resolution order:
    1) dataset_id (if provided and visible to user)
    2) session_data.adata_path
    3) builtin fallback handled in _load_adata_cached when None
    """
    if dataset_id:
        visible = _get_visible_dataset_paths(session_data)
        adata_path = visible.get(dataset_id)
        if not adata_path:
            raise HTTPException(
                status_code=404,
                detail=f"Dataset '{dataset_id}' not found or not visible",
            )
        return adata_path
    return session_data.adata_path


def _extract_gene_min_max(adata: sc.AnnData, gene: str) -> Optional[Dict[str, float]]:
    """
    Return min/max for one gene from an AnnData object.
    Returns None if gene is not found.
    """
    if gene not in adata.var_names:
        return None

    vector = adata[:, gene].X
    if sparse.issparse(vector):
        min_val = float(vector.min())
        max_val = float(vector.max())
    else:
        arr = np.asarray(vector).ravel()
        if arr.size == 0:
            return None
        min_val = float(np.min(arr))
        max_val = float(np.max(arr))

    return {"min": min_val, "max": max_val}


def _get_cached_or_compute_gene_stats(dataset_id: str, adata_path: str, gene: str) -> Optional[Dict[str, float]]:
    """
    Return min/max for (dataset, gene), using TTL cache.
    """
    now = datetime.utcnow()
    cache_key = (dataset_id, _normalize_gene_key(gene))
    cached = GENE_STATS_CACHE.get(cache_key)
    if cached and cached.get("expires_at") and cached["expires_at"] > now:
        return {"min": float(cached["min"]), "max": float(cached["max"])}

    adata = _load_adata_cached(adata_path)
    stats = _extract_gene_min_max(adata, gene)
    if stats is None:
        return None

    GENE_STATS_CACHE[cache_key] = {
        "min": stats["min"],
        "max": stats["max"],
        "expires_at": now + timedelta(seconds=GENE_STATS_CACHE_TTL_SECONDS),
    }
    return stats


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

            # InMemoryBackend doesn't expose _data directly, so we skip aggressive cleanup.
            # Sessions will naturally be garbage collected when no longer referenced.
            # This is a fallback that won't error out.
            pass
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

def build_upload_request(
    email: Optional[str],
    dataset: str,
    saved_files: Dict[str, Optional[str]],
    # Spatial
    spatial_normalization: bool,
    spatial_filtering: bool,
    # Tangram
    use_tangram: bool,
    singlecell_filtering: bool,
    singlecell_normalization: bool,
    gene_selection_mode: Optional[str],
    # Multiome
    use_multiome: bool,
    # Scores
    score_network: bool,
    score_squidpy: bool,
    score_liana_plus: bool,
    score_chromVar: bool,
    score_differential_motif_activity: bool,
    score_motif_enrichment: bool,
    score_FootprintingBias: bool,
    genome: Optional[str],
    # Network Algos
    alg_viper: bool,
    alg_aucell: bool,
    alg_gsva: bool,
    alg_ssgsea: bool,
    # SPONGE-Params
    net_m_score_threshold: Optional[float],
    net_p_adjust: Optional[str],
    net_ensembl_id_col: Optional[str],
    net_feature_col: Optional[str],
    net_rna_types: Optional[str],
    net_max_modules: Optional[int],
    # GENIE3-Params
    genie3_top_n_weights: Optional[int],
    genie3_n_regulatory_genes: Optional[int],
    genie3_n_regulons: Optional[int],
    # Squidpy Flags + Params
    squidpy_moranI: bool,
    squidpy_moranI_n_perms: Optional[int],
    squidpy_moranI_two_tailed: bool,
    squidpy_moranI_corr_method: Optional[str],
    squidpy_gearyC: bool,
    squidpy_gearyC_n_perms: Optional[int],
    squidpy_gearyC_two_tailed: bool,
    squidpy_gearyC_corr_method: Optional[str],
    squidpy_centrality_score: bool,
    squidpy_centrality_score_cluster_key: Optional[str],
    squidpy_co_occurrence: bool,
    squidpy_co_occurrence_cluster_key: Optional[str],
    squidpy_co_occurrence_interval: Optional[int],
    squidpy_co_occurrence_n_splits: Optional[int],
    squidpy_neighborhood_enrichment: bool,
    squidpy_neighborhood_enrichment_cluster_key: Optional[str],
    squidpy_neighborhood_enrichment_library_key: Optional[str],
    squidpy_neighborhood_enrichment_n_perms: Optional[int],
    # ChromVAR Flags + Params
    chromVar_moranI: bool,
    chromVar_moranI_n_perms: Optional[int],
    chromVar_moranI_two_tailed: Optional[str],
    chromVar_moranI_corr_method: Optional[str],
    chromVar_gearyC: bool,
    chromVar_gearyC_n_perms: Optional[int],
    chromVar_gearyC_two_tailed: Optional[str],
    chromVar_gearyC_corr_method: Optional[str],
    chromVar_differential_motif_activity: bool,
    # LIANA
    liana_composition_column: Optional[str],
) -> UploadRequest:
    """
    Convert individual form parameters into a structured UploadRequest model.
    """
    from models import (
        FilesInput,
        SpatialInput,
        TangramInput,
        MultiomeInput,
        ScoresInput,
        NetworkAlgorithms,
        SpongeParams,
        Genie3Params,
        NetworkConfig,
        SquidpyMoranIParams,
        SquidpyGearyCParams,
        SquidpyCentralityParams,
        SquidpyCoOccurrenceParams,
        SquidpyNeighborhoodEnrichmentParams,
        SquidpyConfig,
        ChromVarMoranIParams,
        ChromVarGearyCParams,
        ChromVarConfig,
        LianaConfig,
        UploadRequest,
    )

    files_input = FilesInput(**saved_files)

    spatial_input = SpatialInput(
        normalization=spatial_normalization,
        filtering=spatial_filtering,
    )

    tangram_input = TangramInput(
        use=use_tangram,
        filtering=singlecell_filtering if use_tangram else None,
        normalization=singlecell_normalization if use_tangram else None,
        gene_selection_mode=GeneSelectionMode(gene_selection_mode) if (use_tangram and gene_selection_mode) else None,
    )

    multiome_input = MultiomeInput(use=use_multiome)

    scores_input = ScoresInput(
        network=score_network,
        squidpy=score_squidpy,
        liana_plus=score_liana_plus,
        chromVar=score_chromVar,
        differential_motif_activity=score_differential_motif_activity,
        motif_enrichment=score_motif_enrichment,
        footprinting=score_FootprintingBias,
    )

    network_config = NetworkConfig(
        algorithms=NetworkAlgorithms(
            viper=alg_viper,
            aucell=alg_aucell,
            gsva=alg_gsva,
            ssgsea=alg_ssgsea,
        ),
        sponge_params=SpongeParams(
            m_score_threshold=net_m_score_threshold,
            p_adjust=net_p_adjust,
            ensembl_id_col=net_ensembl_id_col,
            feature_col=net_feature_col,
            rna_types=net_rna_types,
            max_modules=net_max_modules,
        ),
        genie3_params=Genie3Params(
            top_n_weights=genie3_top_n_weights,
            n_regulatory_genes=genie3_n_regulatory_genes,
            n_regulons=genie3_n_regulons,
        ),
    )

    squidpy_config = SquidpyConfig(
        moranI=squidpy_moranI,
        moranI_params=SquidpyMoranIParams(
            n_perms=squidpy_moranI_n_perms,
            two_tailed=squidpy_moranI_two_tailed,
            corr_method=squidpy_moranI_corr_method,
        ),
        gearyC=squidpy_gearyC,
        gearyC_params=SquidpyGearyCParams(
            n_perms=squidpy_gearyC_n_perms,
            two_tailed=squidpy_gearyC_two_tailed,
            corr_method=squidpy_gearyC_corr_method,
        ),
        centrality_score=squidpy_centrality_score,
        centrality_score_params=SquidpyCentralityParams(
            cluster_key=squidpy_centrality_score_cluster_key,
        ),
        co_occurrence=squidpy_co_occurrence,
        co_occurrence_params=SquidpyCoOccurrenceParams(
            cluster_key=squidpy_co_occurrence_cluster_key,
            interval=squidpy_co_occurrence_interval,
            n_splits=squidpy_co_occurrence_n_splits,
        ),
        neighborhood_enrichment=squidpy_neighborhood_enrichment,
        neighborhood_enrichment_params=SquidpyNeighborhoodEnrichmentParams(
            cluster_key=squidpy_neighborhood_enrichment_cluster_key,
            library_key=squidpy_neighborhood_enrichment_library_key,
            n_perms=squidpy_neighborhood_enrichment_n_perms,
        ),
    )

    chromvar_config = ChromVarConfig(
        moranI=chromVar_moranI,
        moranI_params=ChromVarMoranIParams(
            n_perms=chromVar_moranI_n_perms,
            two_tailed=chromVar_moranI_two_tailed,
            corr_method=chromVar_moranI_corr_method,
        ),
        gearyC=chromVar_gearyC,
        gearyC_params=ChromVarGearyCParams(
            n_perms=chromVar_gearyC_n_perms,
            two_tailed=chromVar_gearyC_two_tailed,
            corr_method=chromVar_gearyC_corr_method,
        ),
        differential_motif_activity=chromVar_differential_motif_activity,
    )

    liana_config = LianaConfig(
        composition_column=liana_composition_column,
    )

    return UploadRequest(
        email=email,
        dataset=DatasetType(dataset),
        files=files_input,
        spatial=spatial_input,
        tangram=tangram_input,
        multiome=multiome_input,
        scores=scores_input,
        network=network_config,
        squidpy=squidpy_config,
        chromVar=chromvar_config,
        liana=liana_config,
        genome=genome,
    )


# Main endpoint for receiving form fields and files (multipart/form-data).
# It creates a dedicated job directory, saves all provided files there,
# stores a config JSON for reproducibility, and returns a summary payload.

@app.post("/api/upload", dependencies=[Depends(cookie)])
async def upload(
    # --- metadata ---
    email: Optional[str] = Form(None),
    dataset: str = Form(...),
    # --- files ---
    spatial_h5ad: UploadFile = File(...),
    single_cell_h5ad: Optional[UploadFile] = File(None),
    multiome_rds: Optional[UploadFile] = File(None),
    fragments_tsv_gz: Optional[UploadFile] = File(None),
    fragments_tsv_gz_tbi: Optional[UploadFile] = File(None),
    genie3_network: Optional[UploadFile] = File(None),
    sponge_networkanalysis: Optional[UploadFile] = File(None),
    sponge_networkinteractions: Optional[UploadFile] = File(None),
    liana_genie3_network: Optional[UploadFile] = File(None),
    liana_pathway_network: Optional[UploadFile] = File(None),
    # --- spatial options ---
    spatial_normalization: bool = Form(False),
    spatial_filtering: bool = Form(False),
    # --- tangram options ---
    use_tangram: bool = Form(False),
    singlecell_filtering: bool = Form(False),
    singlecell_normalization: bool = Form(False),
    gene_selection_mode: Optional[str] = Form(None),
    # --- multiome options ---
    use_multiome: bool = Form(False),
    # --- scoring ---
    score_network: bool = Form(False),
    score_squidpy: bool = Form(False),
    score_liana_plus: bool = Form(False),
    score_chromVar: bool = Form(False),
    score_differential_motif_activity: bool = Form(False),
    score_motif_enrichment: bool = Form(False),
    score_FootprintingBias: bool = Form(False),
    genome: Optional[str] = Form(None),
    # --- network algorithms ---
    alg_viper: bool = Form(False),
    alg_aucell: bool = Form(False),
    alg_gsva: bool = Form(False),
    alg_ssgsea: bool = Form(False),
    # --- sponge params ---
    net_m_score_threshold: Optional[float] = Form(None),
    net_p_adjust: Optional[str] = Form(None),
    net_ensembl_id_col: Optional[str] = Form(None),
    net_feature_col: Optional[str] = Form(None),
    net_rna_types: Optional[str] = Form(None),
    net_max_modules: Optional[int] = Form(None),
    # --- genie3 params ---
    genie3_top_n_weights: Optional[int] = Form(None),
    genie3_n_regulatory_genes: Optional[int] = Form(None),
    genie3_n_regulons: Optional[int] = Form(None),
    # --- squidpy params ---
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
    # --- chromvar params ---
    chromVar_moranI: bool = Form(False),
    chromVar_moranI_n_perms: Optional[int] = Form(None),
    chromVar_moranI_two_tailed: Optional[str] = Form(None),
    chromVar_moranI_corr_method: Optional[str] = Form(None),
    chromVar_gearyC: bool = Form(False),
    chromVar_gearyC_n_perms: Optional[int] = Form(None),
    chromVar_gearyC_two_tailed: Optional[str] = Form(None),
    chromVar_gearyC_corr_method: Optional[str] = Form(None),
    chromVar_differential_motif_activity: bool = Form(False),
    # --- liana params ---
    liana_composition_column: Optional[str] = Form(None),
    # --- session ---
    session_data: "SessionData" = Depends(verifier),
):
    """
    Main upload endpoint. Accepts multipart form data with files and options.
    Returns a structured response with processed data paths and configuration.
    """
    print("=== Upload endpoint called ===")
    raw_username = session_data.username
    user_safe = _sanitize_filename(raw_username) or "anon"
    job_id = f"job_{int(time.time() * 1000)}"
    job_dir = BASE_UPLOAD_DIR / f"{job_id}_{user_safe}"
    job_dir.mkdir(parents=True, exist_ok=True)

    # 1) Save all uploaded files to disk
    saved_files_dict = {
        "spatial_h5ad": save_file(spatial_h5ad, job_dir),
        "single_cell_h5ad": save_file(single_cell_h5ad, job_dir),
        "multiome_rds": save_file(multiome_rds, job_dir),
        "fragments_tsv_gz": save_file(fragments_tsv_gz, job_dir),
        "fragments_tsv_gz_tbi": save_file(fragments_tsv_gz_tbi, job_dir),
        "genie3_network": save_file(genie3_network, job_dir),
        "sponge_networkanalysis": save_file(sponge_networkanalysis, job_dir),
        "sponge_networkinteractions": save_file(sponge_networkinteractions, job_dir),
        "liana_genie3_network": save_file(liana_genie3_network, job_dir),
        "liana_pathway_network": save_file(liana_pathway_network, job_dir),
    }


    # 2) Handle RDS to H5AD conversion if needed
    if use_multiome and saved_files_dict.get("multiome_rds"):
        print("Converting RDS to H5AD...")
        rds_path = saved_files_dict["multiome_rds"]
        h5ad_path = re.sub(r"\.rds$", ".h5ad", rds_path)
        log_path = Path(h5ad_path).with_suffix(".log")
        with log_path.open("w") as log_file:
            result = subprocess.run(
                [
                    "Rscript",
                    "../backend/rds_to_h5ad.R",
                    "--rds_path", rds_path,
                    "--assay", "RNA",
                    "--h5ad_path", h5ad_path,
                ],
                stdout=log_file,
                stderr=log_file,
                text=True,
                check=False,
            )
        saved_files_dict["single_cell_h5ad"] = h5ad_path

    # 3) Build structured UploadRequest using helper function
    upload_request = build_upload_request(
        email=email,
        dataset=dataset,
        saved_files=saved_files_dict,
        spatial_normalization=spatial_normalization,
        spatial_filtering=spatial_filtering,
        use_tangram=use_tangram,
        singlecell_filtering=singlecell_filtering,
        singlecell_normalization=singlecell_normalization,
        gene_selection_mode=gene_selection_mode,
        use_multiome=use_multiome,
        score_network=score_network,
        score_squidpy=score_squidpy,
        score_liana_plus=score_liana_plus,
        score_chromVar=score_chromVar,
        score_differential_motif_activity=score_differential_motif_activity,
        score_motif_enrichment=score_motif_enrichment,
        score_FootprintingBias=score_FootprintingBias,
        genome=genome,
        alg_viper=alg_viper,
        alg_aucell=alg_aucell,
        alg_gsva=alg_gsva,
        alg_ssgsea=alg_ssgsea,
        net_m_score_threshold=net_m_score_threshold,
        net_p_adjust=net_p_adjust,
        net_ensembl_id_col=net_ensembl_id_col,
        net_feature_col=net_feature_col,
        net_rna_types=net_rna_types,
        net_max_modules=net_max_modules,
        genie3_top_n_weights=genie3_top_n_weights,
        genie3_n_regulatory_genes=genie3_n_regulatory_genes,
        genie3_n_regulons=genie3_n_regulons,
        squidpy_moranI=squidpy_moranI,
        squidpy_moranI_n_perms=squidpy_moranI_n_perms,
        squidpy_moranI_two_tailed=squidpy_moranI_two_tailed,
        squidpy_moranI_corr_method=squidpy_moranI_corr_method,
        squidpy_gearyC=squidpy_gearyC,
        squidpy_gearyC_n_perms=squidpy_gearyC_n_perms,
        squidpy_gearyC_two_tailed=squidpy_gearyC_two_tailed,
        squidpy_gearyC_corr_method=squidpy_gearyC_corr_method,
        squidpy_centrality_score=squidpy_centrality_score,
        squidpy_centrality_score_cluster_key=squidpy_centrality_score_cluster_key,
        squidpy_co_occurrence=squidpy_co_occurrence,
        squidpy_co_occurrence_cluster_key=squidpy_co_occurrence_cluster_key,
        squidpy_co_occurrence_interval=squidpy_co_occurrence_interval,
        squidpy_co_occurrence_n_splits=squidpy_co_occurrence_n_splits,
        squidpy_neighborhood_enrichment=squidpy_neighborhood_enrichment,
        squidpy_neighborhood_enrichment_cluster_key=squidpy_neighborhood_enrichment_cluster_key,
        squidpy_neighborhood_enrichment_library_key=squidpy_neighborhood_enrichment_library_key,
        squidpy_neighborhood_enrichment_n_perms=squidpy_neighborhood_enrichment_n_perms,
        chromVar_moranI=chromVar_moranI,
        chromVar_moranI_n_perms=chromVar_moranI_n_perms,
        chromVar_moranI_two_tailed=chromVar_moranI_two_tailed,
        chromVar_moranI_corr_method=chromVar_moranI_corr_method,
        chromVar_gearyC=chromVar_gearyC,
        chromVar_gearyC_n_perms=chromVar_gearyC_n_perms,
        chromVar_gearyC_two_tailed=chromVar_gearyC_two_tailed,
        chromVar_gearyC_corr_method=chromVar_gearyC_corr_method,
        chromVar_differential_motif_activity=chromVar_differential_motif_activity,
        liana_composition_column=liana_composition_column,
    )

    # 3b) Convert UploadRequest to Params for OOP dataset structure
    try:
        params = Params.from_upload_request(upload_request)
        print(f"✓ Created Params from UploadRequest")
    except Exception as e:
        print(f"⚠ Could not create Params: {e}; will use raw UploadRequest")
        params = None

    # 4) Convert UploadRequest to dict for pipeline processing
    payload = upload_request.model_dump(exclude_none=False)

    # 5) Run analysis pipeline
    out_dir = await calculate_scores_helper(job_dir, payload)

    print(f"\n=== DEBUG upload() ===")
    print(f"dataset: {dataset}")
    print(f"job_dir: {job_dir}")
    print(f"out_dir: {out_dir}")

    # 6) Identify output files
    adata_path = None
    selected_reason = None
    tangram_adata_path = None
    out_files = {}

    if out_dir is not None:
        out_dir = Path(out_dir)
        # Check for xenium cells file first
        xenium_cells_file = out_dir / "xenium_cells_with_grid_scores.h5ad"
        st_scores_file = out_dir / "st_scores.h5ad"

        if xenium_cells_file.exists() and xenium_cells_file.is_file():
            adata_path = str(xenium_cells_file)
            selected_reason = "xenium_cells_with_grid_scores"
        elif st_scores_file.exists() and st_scores_file.is_file():
            adata_path = str(st_scores_file)
            selected_reason = "st_scores"
        else:
            adata_path = str(saved_files_dict.get("spatial_h5ad"))
            selected_reason = "original_spatial"


        tangram_adata_path = str(out_dir / "tangram_results.h5ad") if (out_dir / "tangram_results.h5ad").exists() else None
        geojson_path = os.path.join(job_dir, "hexagons.geojson")
        if geojson_path:
            out_files["geojson_path"] = f"/api/geojson/{job_id}_{user_safe}"
            print("geojson_path:", geojson_path)
            print("geojson input adata:", adata_path)
            print("geojson data_type:", dataset.lower())

            subprocess.run(
                [
                    "python",
                    "../backend/visium_to_geojson.py",
                    "--adata", adata_path,
                    "--outpath", geojson_path,
                    "--data_type", dataset.lower(),
                ]
            )
            print(f"✓ Generated GeoJSON at {geojson_path}")
        geojson_path = str(out_dir / "hexagons.geojson") if (out_dir / "hexagons.geojson").exists() else None

        # Network files
        if Path(out_dir / "genie_network_filtered_st.csv").exists():
            out_files["genie_network_path"] = str(out_dir / "genie_network_filtered_st.csv")
        if Path(out_dir / "sponge_network_filtered_st.csv").exists():
            out_files["sponge_network_path"] = str(out_dir / "sponge_network_filtered_st.csv")

    if adata_path is None:
        adata_path = saved_files_dict.get("spatial_h5ad")
        selected_reason = "fallback_to_input"

    out_files["adata_path"] = adata_path
    if tangram_adata_path:
        out_files["tangram_adata_path"] = tangram_adata_path

    payload["output_files"] = out_files

    # 7) Create and register Dataset object if Params were successfully created
    if params is not None and adata_path:
        try:
            dataset_obj = DatasetFactory.create_dataset(
                params=params,
                dataset_id=f"{job_id}_{user_safe}",
                user=raw_username,
                adata_path=adata_path,
                tangram_adata_path=tangram_adata_path,
                geojson_path=out_files.get("geojson_path"),
                genie_network_path=out_files.get("genie_network_path"),
                sponge_network_path=out_files.get("sponge_network_path"),
            )
            dataset_registry.register_uploaded_dataset(dataset=dataset_obj)
            print(f"✓ Registered Dataset: {dataset_obj.id}")
        except Exception as e:
            print(f"⚠ Could not register Dataset object: {e}")
            import traceback
            traceback.print_exc()

    # 8) Persist config for reproducibility
    (job_dir / f"{job_id}_config.json").write_text(
        json.dumps(payload, indent=2, default=str), encoding="utf-8"
    )

    # 9) Return clean JSON response
    return payload


@app.get("/api/datasets", dependencies=[Depends(cookie)])
async def get_datasets(session_data: SessionData = Depends(verifier)):
    """
    Return all available datasets for the current user + builtin datasets + shared datasets.
    Filters out datasets with missing files.
    Returns Dataset objects serialized as JSON.
    """
    try:
        # Get all datasets (already serialized as dicts)
        all_datasets = dataset_registry.get_all_datasets(as_dict=True)
        print(f"Fetching datasets for user: {session_data.username}")

        # Build response: filter out datasets with missing files
        datasets_json = {
            "builtin": {},
            "uploaded": {},
        }

        # Process builtin datasets
        for dataset_id, dataset_dict in all_datasets.get("builtin", {}).items():
            adata_path = dataset_dict.get("adata_path")
            if adata_path and Path(adata_path).exists():
                datasets_json["builtin"][dataset_id] = dataset_dict
                print(f"✓ Included builtin dataset: {dataset_id}")
            else:
                print(f"✗ Skipped builtin dataset {dataset_id}: missing adata at {adata_path}")

        # Process user's uploaded datasets
        user = session_data.username
        for dataset_id, dataset_dict in all_datasets.get("uploaded", {}).items():
            if dataset_dict.get("user") == user:
                adata_path = dataset_dict.get("adata_path")
                if adata_path and Path(adata_path).exists():
                    datasets_json["uploaded"][dataset_id] = dataset_dict
                    print(f"✓ Included uploaded dataset: {dataset_id}")
                else:
                    print(f"✗ Skipped uploaded dataset {dataset_id}: missing adata at {adata_path}")

        print(f"Returning datasets: {list(datasets_json['builtin'].keys())} builtin + {list(datasets_json['uploaded'].keys())} uploaded")
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

    if genie_network_path.path and not os.path.exists(genie_network_path.path):
        raise HTTPException(status_code=400, detail="Genie network file not found")
    elif genie_network_path.path:
        session_data.genie_network_path = genie_network_path.path
    elif not genie_network_path.path:
        session_data.genie_network_path = None
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

    if sponge_network_path.path and not os.path.exists(sponge_network_path.path):
        raise HTTPException(status_code=400, detail="Sponge network file not found")
    elif sponge_network_path.path:
        session_data.sponge_network_path = sponge_network_path.path
    elif not sponge_network_path.path:
        session_data.sponge_network_path = None

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
async def get_obsm_column(
    table: str,
    column: str,
    dataset_id: Optional[str] = None,
    session_data: SessionData = Depends(verifier),
):
    print(f"[DEBUG] Endpoint called: table={table}, column={column}")
    print(f"[DEBUG] Session username: {session_data.username}")
    adata_path = _resolve_adata_path(session_data, dataset_id)
    adata = _load_adata_cached(adata_path)

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
    barcode: str,
    dataset_id: Optional[str] = None,
    session_data: SessionData = Depends(verifier)
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
    adata_path = _resolve_adata_path(session_data, dataset_id)
    adata = _load_adata_cached(adata_path)
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
    gene: str,
    dataset_id: Optional[str] = None,
    session_data: SessionData = Depends(verifier)
):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:3000/X/ENSG00000241860`
    """
    adata_path = _resolve_adata_path(session_data, dataset_id)
    adata = _load_adata_cached(adata_path)
    expressions = adata[:, gene].X.toarray().flatten().tolist()
    barcodes = adata.obs.index
    return {
        barcode: expression
        for barcode, expression in zip(barcodes, expressions)
    }


@app.get("/X_stats/{gene}", dependencies=[Depends(cookie)])
async def get_X_stats_by_gene(
    gene: str,
    dataset_ids: Optional[str] = None,
    session_data: SessionData = Depends(verifier),
):
    """
    Return aggregated expression min/max for `gene` across selected visible datasets.

    Query params:
    - dataset_ids: optional comma-separated dataset ids.
      If omitted, all visible datasets for the current user are used.

    Response:
    {
      "gene": "...",
      "dataset_ids": [...],
      "global_min": float,
      "global_max": float,
      "per_dataset": {"dataset_id": {"min": ..., "max": ...}},
      "missing_in": [...]
    }
    """
    requested_gene = (gene or "").strip()
    if not requested_gene:
        raise HTTPException(status_code=400, detail="Gene must not be empty")

    visible_paths = _get_visible_dataset_paths(session_data)
    if not visible_paths:
        raise HTTPException(status_code=404, detail="No visible datasets available")

    selected_ids: List[str]
    if dataset_ids:
        raw_ids = [x.strip() for x in dataset_ids.split(",") if x.strip()]
        selected_ids = [dataset_id for dataset_id in raw_ids if dataset_id in visible_paths]
        if not selected_ids:
            raise HTTPException(
                status_code=404,
                detail="None of the requested dataset_ids are visible/available",
            )
    else:
        selected_ids = list(visible_paths.keys())

    per_dataset: Dict[str, Dict[str, float]] = {}
    missing_in: List[str] = []

    for dataset_id in selected_ids:
        adata_path = visible_paths[dataset_id]
        try:
            stats = _get_cached_or_compute_gene_stats(dataset_id, adata_path, requested_gene)
        except Exception as e:
            print(f"[X_stats] Failed for dataset={dataset_id}, gene={requested_gene}: {e}")
            stats = None

        if stats is None:
            missing_in.append(dataset_id)
        else:
            per_dataset[dataset_id] = stats

    if not per_dataset:
        raise HTTPException(
            status_code=404,
            detail=f"Gene '{requested_gene}' not found in any selected dataset",
        )

    global_min = float(min(v["min"] for v in per_dataset.values()))
    global_max = float(max(v["max"] for v in per_dataset.values()))

    return {
        "gene": requested_gene,
        "dataset_ids": selected_ids,
        "global_min": global_min,
        "global_max": global_max,
        "per_dataset": per_dataset,
        "missing_in": missing_in,
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
    print("[DEBUG] get available motifs")
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


@app.get("/api/cell_types", dependencies=[Depends(cookie)])
async def get_available_cell_types(
    session_data: SessionData = Depends(verifier),
    dataset_id: Optional[str] = None,
):
    """
    Return the unique cell types from spot_obj_footprints.rds for the current dataset.
    Uses a small inline Rscript so no extra dependency is needed.
    """
    print("[DEBUG] get available cell types")
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
    print(f"[DEBUG] adata_path = {adata_path}")
    try:
        adata = _load_adata_cached(adata_path)
        cell_types = list(adata.obs["cell_type"].unique())
        print(f"[get_available_cell_types] Found cell types: {cell_types}, dataset_id={dataset_id}")

        return {"cell_types": cell_types}
    except subprocess.TimeoutExpired:
        raise HTTPException(status_code=504, detail="Rscript timed out reading cell types")


@app.post("/api/compute_footprint", dependencies=[Depends(cookie)])
async def compute_footprint(
    motif: List[str] = Form(...),
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
    motifs_csv = ",".join(motif)
    cmd = [
        "Rscript", "../backend/calc_multiome_scores/compute_additional_footprints.R",
        "--outdir", str(out_dir),
        "--motifs", motifs_csv,
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

    # Collect one result entry per requested motif
    results = []
    all_datasets = dataset_registry.get_all_datasets()
    for single_motif in motif:
        pdf_filename = f"footprint_{single_motif}_{cluster_by}.pdf"
        pdf_path = out_dir / pdf_filename
        if not pdf_path.exists():
            raise HTTPException(status_code=500, detail=f"Expected output {pdf_filename} was not created")

        relative_pdf = str(pdf_path.relative_to(BASE_UPLOAD_DIR.resolve()))
        for category in ("uploaded", "builtin"):
            for ds_id, info in all_datasets.get(category, {}).items():
                if info.get("adata_path") == adata_path:
                    existing = info.get("footprint_list") or []
                    if relative_pdf not in existing:
                        dataset_registry.update_dataset_paths(ds_id, footprint_list=existing + [relative_pdf])
                    break

        results.append({
            "footprint_url": f"/api/download/{relative_pdf}",
            "relative_path": relative_pdf,
            "motif": single_motif,
            "cluster_by": cluster_by,
        })

    return {"results": results}


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=3000)
    # for merit
    # uvicorn.run(app, host="0.0.0.0", port=3005)
