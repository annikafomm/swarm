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
from enum import Enum
from pathlib import Path
from typing import Any, Dict, Optional
from uuid import UUID, uuid4

import numpy as np
import pandas as pd
import scanpy as sc
import uvicorn
from app import calculate_scores_helper

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
from fastapi.responses import FileResponse
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
BASE_UPLOAD_DIR = Path.cwd() / "../backend/uploads"
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

class GeneSelectionMode(str, Enum):
    ctg = "ctg"
    hvg = "hvg"
    spapros = "spapros"
    svg = "svg"
    none = "None"


class Method(str, Enum):
    Genie3 = "Genie3"
    Sponge = "Sponge"

class BaseModel(PydanticBaseModel):
    class Config:
        arbitrary_types_allowed = True


class SessionData(BaseModel):
    username: str
    adata_path: str = None
    adata: sc.AnnData | None = None
    genie_network_path: str | None = None
    sponge_network_path: str | None = None



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
    # Scores
    score_network: bool = Form(False),
    score_squidpy: bool = Form(False),
    score_liana_plus: bool = Form(False),
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
    # LIANA
    liana_composition_column: Optional[str] = Form(None),
    liana_genie3_network: Optional[UploadFile] = File(None),
    liana_pathway_network: Optional[UploadFile] = File(None),
    session_data: "SessionData" = Depends(verifier),
):

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

    raw_username = session_data.username
    user_safe = _sanitize_filename(raw_username) or "anon"

    job_id = f"job_{int(time.time() * 1000)}"
    job_dir = BASE_UPLOAD_DIR / f"{job_id}_{user_safe}"
    job_dir.mkdir(parents=True, exist_ok=True)

    # 3) Save all files into that directory
    # Dateien speichern (alles in ./uploads)
    saved_files = {
        "spatial_h5ad": save_file(spatial_h5ad, job_dir),
        "single_cell_h5ad": save_file(single_cell_h5ad, job_dir),
        "genie3_network": save_file(genie3_network, job_dir),
        "sponge_networkanalysis": save_file(sponge_networkanalysis, job_dir),
        "sponge_networkinteractions": save_file(
            sponge_networkinteractions, job_dir
        ),
        "liana_genie3_network": save_file(liana_genie3_network, job_dir),
        "liana_pathway_network": save_file(liana_pathway_network, job_dir),
    }

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
        "scores": {
            "network": score_network,
            "squidpy": score_squidpy,
            "liana_plus": score_liana_plus,
        },
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

        print("selected adata_path:", adata_path)
        print("selection reason:", selected_reason)




    if adata_path is None:
        adata_path = saved_files.get("spatial_h5ad")
        selected_reason = "Fallback to uploaded spatial_h5ad (no pipeline output found)"
        print("WARNING:", selected_reason)

    if out_dir is None:
        out_dir = job_dir

    if adata_path is not None:
        out_files["adataPath"] = adata_path

        data_type = dataset.value.lower()

        geojson_out = os.path.join(out_dir, "hexagons.geojson")
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
                os.path.join(out_dir, "hexagons.geojson"),
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


@app.get("/api/hexagon/{user}/{subdir}/{filename}")
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

    session_data.adata_path = adata_path.path
    session_data.adata = adata
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

    # Get the geneset from the name
    gene_set = session_data.adata.uns["genie_genesets"].get(
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
    gene_set = session_data.adata.uns["sponge_genesets"].get(
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
    available_scores = session_data.adata.obsm.keys()
    available_scores = [score for score in available_scores if score.endswith("_genie3") or score.endswith("_sponge")]

    row_data = {}
    for score in available_scores:
        obsm_data = session_data.adata.obsm[score]
        if not isinstance(obsm_data, pd.DataFrame):
            obsm_data = pd.DataFrame(
                obsm_data, index=session_data.adata.obs_names
            )
        row_data[score] = obsm_data.loc[barcode].to_dict()
    return row_data


@app.get("/X/{gene}", dependencies=[Depends(cookie)])
async def get_X_by_gene(
    gene: str, session_data: SessionData = Depends(verifier)
):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:3000/X/ENSG00000241860`
    """
    expressions = session_data.adata[:, gene].X.toarray().flatten().tolist()
    barcodes = session_data.adata.obs.index
    return {
        barcode: expression
        for barcode, expression in zip(barcodes, expressions)
    }


if __name__ == "__main__":
    uvicorn.run(app, host="0.0.0.0", port=3000)
