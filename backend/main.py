from uuid import UUID, uuid4

import pandas as pd
import scanpy as sc
import uvicorn
from fastapi import Depends, FastAPI, HTTPException, Response
from fastapi.middleware.cors import CORSMiddleware
from fastapi_sessions.backends.implementations import InMemoryBackend
from fastapi_sessions.frontends.implementations import (
    CookieParameters,
    SessionCookie,
)
from fastapi_sessions.session_verifier import SessionVerifier
from pydantic import BaseModel as PydanticBaseModel


class BaseModel(PydanticBaseModel):
    class Config:
        arbitrary_types_allowed = True


class SessionData(BaseModel):
    username: str
    adata: sc.AnnData | None = None


cookie_params = CookieParameters(secure=False, httponly=True, samesite="lax")

# Uses UUID
cookie = SessionCookie(
    cookie_name="cookie",
    identifier="general_verifier",
    auto_error=True,
    secret_key="DONOTUSE",
    cookie_params=cookie_params,
)
backend = InMemoryBackend[UUID, SessionData]()


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


verifier = BasicVerifier(
    identifier="general_verifier",
    auto_error=True,
    backend=backend,
    auth_http_exception=HTTPException(
        status_code=403, detail="invalid session"
    ),
)

app = FastAPI()

origins = [
    "http://localhost",
    "http://127.0.0.1",
    "http://localhost:4200",
    "http://127.0.0.1:4200",
]


app.add_middleware(
    CORSMiddleware,
    allow_origins=origins,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)


@app.get("/")
async def read_root():
    return {"message": "Hello, FastAPI!"}


@app.post("/create_session/{name}")
async def create_session(name: str, response: Response):
    """
    Example: `curl -c cookies.txt -X POST http://127.0.0.1:8000/create_session/mopitas`
    """
    session = uuid4()
    data = SessionData(username=name)

    await backend.create(session, data)
    cookie.attach_to_response(response, session)

    return f"created session for {name}"


@app.get("/whoami", dependencies=[Depends(cookie)])
async def whoami(session_data: SessionData = Depends(verifier)):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:8000/whoami`
    """
    return session_data.username


@app.post("/delete_session")
async def del_session(response: Response, session_id: UUID = Depends(cookie)):
    """
    Example: `curl -b cookies.txt -X POST http://127.0.0.1:8000/delete_session`
    """
    await backend.delete(session_id)
    cookie.delete_from_response(response)
    return "deleted session"


class AnnDataPath(BaseModel):
    path: str


@app.post("/read_adata")
async def read_adata(
    adata_path: AnnDataPath, session_id: UUID = Depends(cookie)
):
    """
    Example:
    ```
    curl -b cookies.txt \
    -H "Content-Type: application/json" \
    -X POST http://127.0.0.1:8000/read_adata \
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
    Example: `curl -b cookies.txt http://127.0.0.1:8000/obs/cell_type`
    """
    return session_data.adata.obs[column].to_dict()


@app.get("/var/{column}", dependencies=[Depends(cookie)])
async def get_var_column(
    column: str, session_data: SessionData = Depends(verifier)
):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:8000/var/n_cells`
    """
    return session_data.adata.var[column].to_dict()


@app.get("/obsm/{table}/{column}", dependencies=[Depends(cookie)])
async def get_obsm_column(
    table: str, column: str, session_data: SessionData = Depends(verifier)
):
    """
    Example: `curl -b cookies.txt http://127.0.0.1:8000/obsm/ligand_receptor_cosine_similarity/LGALS9^PTPRC`
    """
    return session_data.adata.obsm[table][column].to_dict()


if __name__ == "__main__":
    app.state.data = None
    uvicorn.run(app, host="0.0.0.0", port=8000)
