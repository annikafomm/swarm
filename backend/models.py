"""
Centralized Pydantic models for API requests and responses.
Replaces scattered Form() parameters with hierarchical, validated data structures.
"""

from enum import Enum
from typing import Any, Dict, List, Optional
from datetime import datetime
from pydantic import BaseModel, Field, field_validator


# =============================================================================
# Enumerations
# =============================================================================

class DatasetType(str, Enum):
    """Supported spatial dataset types."""
    Visium = "Visium"
    Xenium = "Xenium"
    Multiome = "Multiome"


class GeneSelectionMode(str, Enum):
    """Gene selection strategies for preprocessing."""
    ctg = "ctg"
    hvg = "hvg"
    spapros = "spapros"
    svg = "svg"
    none = "None"


class Genome(str, Enum):
    """Reference genome versions."""
    hg37 = "hg37"
    hg38 = "hg38"


class Method(str, Enum):
    """Network inference methods."""
    Genie3 = "Genie3"
    Sponge = "Sponge"


# =============================================================================
# Base Classes
# =============================================================================

class BaseConfig(BaseModel):
    """Base configuration with common settings."""
    class Config:
        arbitrary_types_allowed = True
        use_enum_values = False  # Keep enum objects, not strings


# =============================================================================
# File Inputs
# =============================================================================

class FilesInput(BaseConfig):
    """Consolidated file upload paths (after saving to disk)."""
    spatial_h5ad: str
    single_cell_h5ad: Optional[str] = None
    multiome_rds: Optional[str] = None
    fragments_tsv_gz: Optional[str] = None
    fragments_tsv_gz_tbi: Optional[str] = None
    genie3_network: Optional[str] = None
    sponge_networkanalysis: Optional[str] = None
    sponge_networkinteractions: Optional[str] = None
    liana_genie3_network: Optional[str] = None
    liana_pathway_network: Optional[str] = None


# =============================================================================
# Preprocessing Inputs
# =============================================================================

class SpatialInput(BaseConfig):
    """Spatial data processing options."""
    normalization: bool = False
    filtering: bool = False


class TangramInput(BaseConfig):
    """Single-cell to spatial mapping (Tangram) configuration."""
    use: bool = False
    filtering: Optional[bool] = None
    normalization: Optional[bool] = None
    gene_selection_mode: Optional[GeneSelectionMode] = None

    @field_validator("filtering", "normalization", mode="before")
    @classmethod
    def only_if_used(cls, v, info):
        """Parameters only valid if use=True."""
        if not info.data.get("use") and v is not None:
            return None
        return v


class MultiomeInput(BaseConfig):
    """Multiome (ATAC-seq) specific configuration."""
    use: bool = False


# =============================================================================
# Scoring Configurations
# =============================================================================

class ScoresInput(BaseConfig):
    """Overall scoring enablement flags."""
    network: bool = False
    squidpy: bool = False
    liana_plus: bool = False
    chromVar: bool = False
    differential_motif_activity: bool = False
    motif_enrichment: bool = False
    footprinting: bool = False


class NetworkAlgorithms(BaseConfig):
    """Network inference algorithms."""
    viper: bool = False
    aucell: bool = False
    gsva: bool = False
    ssgsea: bool = False


class SpongeParams(BaseConfig):
    """SPONGE network analysis parameters."""
    m_score_threshold: Optional[float] = None
    p_adjust: Optional[str] = None
    ensembl_id_col: Optional[str] = None
    feature_col: Optional[str] = None
    rna_types: Optional[str] = None
    max_modules: Optional[int] = None


class Genie3Params(BaseConfig):
    """GENIE3 network inference parameters."""
    top_n_weights: Optional[int] = None
    n_regulatory_genes: Optional[int] = None
    n_regulons: Optional[int] = None


class NetworkConfig(BaseConfig):
    """Complete network analysis configuration."""
    algorithms: NetworkAlgorithms = Field(default_factory=NetworkAlgorithms)
    sponge_params: SpongeParams = Field(default_factory=SpongeParams)
    genie3_params: Genie3Params = Field(default_factory=Genie3Params)


# =============================================================================
# Spatial Analysis Algorithms
# =============================================================================

class SquidpyMoranIParams(BaseConfig):
    """Moran's I autocorrelation test parameters."""
    n_perms: Optional[int] = None
    two_tailed: bool = False
    corr_method: Optional[str] = None


class SquidpyGearyCParams(BaseConfig):
    """Geary's C autocorrelation test parameters."""
    n_perms: Optional[int] = None
    two_tailed: bool = False
    corr_method: Optional[str] = None


class SquidpyCentralityParams(BaseConfig):
    """Centrality score parameters."""
    cluster_key: Optional[str] = None


class SquidpyCoOccurrenceParams(BaseConfig):
    """Co-occurrence analysis parameters."""
    cluster_key: Optional[str] = None
    interval: Optional[int] = None
    n_splits: Optional[int] = None


class SquidpyNeighborhoodEnrichmentParams(BaseConfig):
    """Neighborhood enrichment parameters."""
    cluster_key: Optional[str] = None
    library_key: Optional[str] = None
    n_perms: Optional[int] = None


class SquidpyConfig(BaseConfig):
    """Complete Squidpy spatial analysis configuration."""
    moranI: bool = False
    moranI_params: SquidpyMoranIParams = Field(default_factory=SquidpyMoranIParams)
    gearyC: bool = False
    gearyC_params: SquidpyGearyCParams = Field(default_factory=SquidpyGearyCParams)
    centrality_score: bool = False
    centrality_score_params: SquidpyCentralityParams = Field(default_factory=SquidpyCentralityParams)
    co_occurrence: bool = False
    co_occurrence_params: SquidpyCoOccurrenceParams = Field(default_factory=SquidpyCoOccurrenceParams)
    neighborhood_enrichment: bool = False
    neighborhood_enrichment_params: SquidpyNeighborhoodEnrichmentParams = Field(default_factory=SquidpyNeighborhoodEnrichmentParams)


# =============================================================================
# ChromVAR Configuration
# =============================================================================

class ChromVarMoranIParams(BaseConfig):
    """ChromVAR Moran's I parameters."""
    n_perms: Optional[int] = None
    two_tailed: Optional[str] = None  # "oneTailed" or "twoTailed"
    corr_method: Optional[str] = None


class ChromVarGearyCParams(BaseConfig):
    """ChromVAR Geary's C parameters."""
    n_perms: Optional[int] = None
    two_tailed: Optional[str] = None
    corr_method: Optional[str] = None


class ChromVarConfig(BaseConfig):
    """Complete ChromVAR configuration."""
    moranI: bool = False
    moranI_params: ChromVarMoranIParams = Field(default_factory=ChromVarMoranIParams)
    gearyC: bool = False
    gearyC_params: ChromVarGearyCParams = Field(default_factory=ChromVarGearyCParams)
    differential_motif_activity: bool = False


# =============================================================================
# LIANA Configuration
# =============================================================================

class LianaConfig(BaseConfig):
    """LIANA ligand-receptor analysis configuration."""
    composition_column: Optional[str] = None


# =============================================================================
# Request/Response Models
# =============================================================================

class UploadRequest(BaseConfig):
    """
    Complete upload request configuration.
    Replaces 90+ scattered Form() parameters with hierarchical structure.
    """
    # Core metadata
    email: Optional[str] = None
    dataset: DatasetType = Field(..., description="Dataset type (Visium or Xenium)")

    # Files (populated after saving uploads to disk)
    files: FilesInput = Field(..., description="File paths after saving")

    # Processing options
    spatial: SpatialInput = Field(default_factory=SpatialInput)
    tangram: TangramInput = Field(default_factory=TangramInput)
    multiome: MultiomeInput = Field(default_factory=MultiomeInput)

    # Scores and analysis
    scores: ScoresInput = Field(default_factory=ScoresInput)
    network: NetworkConfig = Field(default_factory=NetworkConfig)
    squidpy: SquidpyConfig = Field(default_factory=SquidpyConfig)
    chromVar: ChromVarConfig = Field(default_factory=ChromVarConfig)
    liana: LianaConfig = Field(default_factory=LianaConfig)

    # Genome reference
    genome: Optional[str] = None


class OutputFiles(BaseConfig):
    """Output files generated by analysis pipeline."""
    adata_path: Optional[str] = None
    tangram_adata_path: Optional[str] = None
    geojson_path: Optional[str] = None
    genie_network_path: Optional[str] = None
    sponge_network_path: Optional[str] = None
    footprint_list: Optional[List[str]] = None


class UploadResponse(BaseConfig):
    """Standardized upload response structure."""
    email: Optional[str] = None
    dataset: str
    spatial: Dict[str, Any]
    files: Dict[str, Any]
    tangram: Dict[str, Any]
    multiome: Dict[str, Any]
    scores: Dict[str, Any]
    genome: Optional[str] = None
    network: Dict[str, Any]
    squidpy: Dict[str, Any]
    chromVar: Dict[str, Any]
    liana: Dict[str, Any]
    output_files: Optional[OutputFiles] = None
    created_at: datetime = Field(default_factory=datetime.now)
