"""
Dataset structure module providing OOP abstraction for datasets and parameters.

Provides:
- Params: Parsed and validated parameters from upload request
- Dataset: Abstract base class for all dataset types
- VisiumDataset: Concrete implementation for Visium spatial transcriptomics
- DatasetFactory: Factory for creating appropriate dataset types
"""

from __future__ import annotations
from abc import ABC, abstractmethod
from pathlib import Path
from typing import Optional, Dict, Any
from datetime import datetime
from enum import Enum

from models import UploadRequest, DatasetType


# =============================================================================
# Params Class - Validated Parameter Container
# =============================================================================

class Params(UploadRequest):
    """
    Extended UploadRequest with validation and metadata extraction methods.
    Represents all parameters needed to create and configure a dataset.
    """

    @classmethod
    def from_upload_request(cls, upload_request: UploadRequest) -> Params:
        """
        Convert UploadRequest to Params (with validation).

        Args:
            upload_request: Raw UploadRequest from form submission

        Returns:
            Params: Validated and parsed parameters

        Raises:
            ValueError: If parameters are invalid or missing
        """
        # Since Params extends UploadRequest, we can just convert via model_dump
        data = upload_request.model_dump()
        params = cls(**data)
        params.validate_file_paths()
        return params

    def validate_file_paths(self) -> None:
        """
        Validate that uploaded file paths exist and are accessible.

        Raises:
            ValueError: If critical files are missing or inaccessible
        """
        if self.files.spatial_h5ad:
            path = Path(self.files.spatial_h5ad)
            if not path.exists():
                raise ValueError(f"Spatial H5AD file not found: {path}")
            if not path.is_file():
                raise ValueError(f"Spatial H5AD is not a file: {path}")

    def get_dataset_type(self) -> str:
        """Extract dataset type string (for pipeline use)."""
        return self.dataset.value if isinstance(self.dataset, Enum) else str(self.dataset)

    def get_registry_dataset_type(self) -> str:
        """
        Get the dataset type for registry/storage purposes.
        Returns 'multiome' if multiome analysis is enabled, otherwise the spatial type.
        """
        if self.is_multiome_enabled():
            return "multiome"
        return self.get_dataset_type()
    def get_dataset_alias(self) -> str:
        """Generate user-friendly dataset alias."""
        return f"{self.get_dataset_type()} Dataset"

    def is_tangram_enabled(self) -> bool:
        """Check if tangram analysis is enabled."""
        return self.tangram.use if self.tangram else False

    def is_multiome_enabled(self) -> bool:
        """Check if multiome (ATAC-seq) analysis is enabled."""
        return self.multiome.use if self.multiome else False

    def has_network_scores(self) -> bool:
        """Check if any network scoring is enabled."""
        return (
            self.network.algorithms.viper
            or self.network.algorithms.aucell
            or self.network.algorithms.gsva
            or self.network.algorithms.ssgsea
        ) if self.network else False


# =============================================================================
# Dataset Classes - OOP Dataset Abstraction
# =============================================================================

class Dataset(ABC):
    """
    Abstract base class for all dataset types.
    Encapsulates dataset metadata and paths.
    """

    def __init__(
        self,
        dataset_id: str,
        alias: str,
        adata_path: str,
        user: str = "anonymous",
        params: Optional[Params] = None,
        tangram_adata_path: Optional[str] = None,
        geojson_path: Optional[str] = None,
        genie_network_path: Optional[str] = None,
        sponge_network_path: Optional[str] = None,
        tf_graph_path: Optional[str] = None,
        created_at: Optional[datetime] = None,
        **metadata
    ):
        """
        Initialize a Dataset.

        Args:
            dataset_id: Unique identifier (e.g., "job_123_user")
            alias: User-friendly name
            adata_path: Path to main AnnData H5AD file
            user: User who created/owns this dataset
            params: Original Params object that created this dataset
            tangram_adata_path: Optional path to tangram-processed data
            geojson_path: Optional path or URL to GeoJSON visualization
            genie_network_path: Optional path to GENIE3 network
            sponge_network_path: Optional path to SPONGE network
            tf_graph_path: Optional path to TF graph JSON file
            created_at: Creation timestamp
            **metadata: Additional dataset-specific metadata
        """
        self.id = dataset_id
        self.alias = alias
        self._adata_path = adata_path
        self.user = user
        self.params = params
        self.tangram_adata_path = tangram_adata_path
        self.geojson_path = geojson_path
        self.genie_network_path = genie_network_path
        self.sponge_network_path = sponge_network_path
        self.tf_graph_path = tf_graph_path
        self.created_at = created_at or datetime.utcnow()
        self.metadata = metadata

    @property
    def type(self) -> str:
        """Dataset type (e.g., 'visium', 'xenium')"""
        return self.__class__.__name__.replace("Dataset", "").lower()

    @property
    def adata_path(self) -> str:
        """Primary AnnData file path."""
        return self._adata_path

    def validate_paths(self) -> bool:
        """
        Verify that all referenced files exist.

        Returns:
            True if all required paths are valid, False otherwise
        """
        required_path = Path(self.adata_path)
        if not required_path.exists():
            print(f"✗ Dataset {self.id}: adata not found at {self.adata_path}")
            return False

        optional_paths = [self.tangram_adata_path, self.genie_network_path, self.sponge_network_path]
        for path_str in optional_paths:
            if path_str:
                path = Path(path_str)
                if not path.exists():
                    print(f"⚠ Dataset {self.id}: optional file missing at {path_str}")

        return True

    def to_dict(self) -> Dict[str, Any]:
        """
        Serialize dataset to dictionary (for JSON responses).

        Returns:
            Dictionary representation suitable for JSON serialization
        """
        return {
            "id": self.id,
            "type": self.type,
            "alias": self.alias,
            "user": self.user,
            "adata_path": self.adata_path,
            "tangram_adata_path": self.tangram_adata_path,
            "geojson_path": self.geojson_path,
            "genie_network_path": self.genie_network_path,
            "sponge_network_path": self.sponge_network_path,
            "tf_graph_path": self.tf_graph_path,
            "created_at": self.created_at.isoformat() if isinstance(self.created_at, datetime) else self.created_at,
            "metadata": self.metadata,
        }

    def to_registry_dict(self) -> Dict[str, Any]:
        """
        Serialize to registry-compatible format (for JSON persistence).
        Used by DatasetRegistry to store in JSON.

        Returns:
            Dictionary compatible with registry JSON format
        """
        return self.to_dict()

    @classmethod
    def from_registry_dict(cls, data: Dict[str, Any]) -> Dataset:
        """
        Deserialize from registry dictionary.
        Routes to appropriate subclass based on 'type' field using DatasetFactory.

        Args:
            data: Dictionary from registry (should have 'type' field)

        Returns:
            Appropriate Dataset subclass instance (VisiumDataset, XeniumDataset, or MultiomeDataset)
        """
        # Delegate to DatasetFactory which properly routes by type field
        # This ensures multiome and xenium datasets are deserialized to correct classes
        from dataset_structure import DatasetFactory
        return DatasetFactory.from_registry_dict(data)

    @classmethod
    @abstractmethod
    def from_params(cls, params: Params, dataset_id: str, user: str = "anonymous") -> Dataset:
        """
        Factory method: Create Dataset from Params after upload/pipeline completion.

        Args:
            params: Validated parameters from upload
            dataset_id: Generated dataset ID (e.g., "job_123_user")
            user: User who owns this dataset

        Returns:
            Dataset instance of appropriate type
        """
        pass

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(id={self.id}, type={self.type}, alias={self.alias})"


# =============================================================================
# VisiumDataset - Visium-Specific Implementation
# =============================================================================

class VisiumDataset(Dataset):
    """
    Concrete implementation for Visium spatial transcriptomics datasets.
    """

    def __init__(
        self,
        dataset_id: str,
        alias: str,
        adata_path: str,
        user: str = "anonymous",
        params: Optional[Params] = None,
        tangram_adata_path: Optional[str] = None,
        geojson_path: Optional[str] = None,
        genie_network_path: Optional[str] = None,
        sponge_network_path: Optional[str] = None,
        tf_graph_path: Optional[str] = None,
        created_at: Optional[datetime] = None,
        dataset_type: Optional[str] = None,  # e.g., "Visium", "Xenium"
        use_tangram: bool = False,
        use_multiome: bool = False,
        adata_st_scores_path: Optional[str] = None,  # Original ST data with scores
        adata_tg_scores_path: Optional[str] = None,  # Tangram-projected data with scores
        **metadata
    ):
        """
        Initialize a VisiumDataset.

        Args:
            dataset_type: Original dataset type from params (e.g., "Visium", "XeniumRIF1")
            use_tangram: Whether tangram deconvolution was applied
            use_multiome: Whether multiome data was processed
            adata_st_scores_path: Path to original ST data with computed scores (optional, when scores available)
            adata_tg_scores_path: Path to tangram-projected data with scores (optional, when tangram+scores)
            tf_graph_path: Optional path to TF graph JSON file
            **metadata: Additional Visium-specific metadata
        """
        super().__init__(
            dataset_id=dataset_id,
            alias=alias,
            adata_path=adata_path,
            user=user,
            params=params,
            tangram_adata_path=tangram_adata_path,
            geojson_path=geojson_path,
            genie_network_path=genie_network_path,
            sponge_network_path=sponge_network_path,
            tf_graph_path=tf_graph_path,
            created_at=created_at,
            **metadata
        )
        self.dataset_type = dataset_type  # e.g., "Visium"
        self.use_tangram = use_tangram
        self.use_multiome = use_multiome
        self.adata_st_scores_path = adata_st_scores_path  # For datasets with computed scores
        self.adata_tg_scores_path = adata_tg_scores_path  # For tangram+scores datasets

    def to_dict(self) -> Dict[str, Any]:
        """Serialize VisiumDataset to dictionary."""
        data = super().to_dict()
        data.update({
            "dataset_type": self.dataset_type,
            "use_tangram": self.use_tangram,
            "use_multiome": self.use_multiome,
            "adata_st_scores_path": self.adata_st_scores_path,
            "adata_tg_scores_path": self.adata_tg_scores_path,
            "tf_graph_path": self.tf_graph_path,
        })
        return data

    @classmethod
    def from_params(
        cls,
        params: Params,
        dataset_id: str,
        user: str = "anonymous",
        adata_path: Optional[str] = None,
        tangram_adata_path: Optional[str] = None,
        geojson_path: Optional[str] = None,
        genie_network_path: Optional[str] = None,
        sponge_network_path: Optional[str] = None,
        # Xenium-specific (accepted but ignored)
        xenium_grid_adata_path: Optional[str] = None,
        # Multiome-specific (accepted but ignored)
        adata_st_scores_path: Optional[str] = None,
        adata_tg_scores_path: Optional[str] = None,
        adata_map_path: Optional[str] = None,
        adata_map_X_csv_path: Optional[str] = None,
        adata_map_var_csv_path: Optional[str] = None,
        calc_scores_log_path: Optional[str] = None,
        global_motif_analysis_path: Optional[str] = None,
        motif_to_tf_csv_path: Optional[str] = None,
        spot_obj_chromvar_path: Optional[str] = None,
        spot_obj_footprints_path: Optional[str] = None,
        dissociated_obj_footprints_path: Optional[str] = None,
        chromvar_scores_csv_path: Optional[str] = None,
        diff_motif_activity_csv_paths: Optional[Dict[str, str]] = None,
        footprint_pdf_paths: Optional[Dict[str, str]] = None,
    ) -> VisiumDataset:
        """
        Create VisiumDataset from Params and output paths.

        Args:
            params: Validated parameters
            dataset_id: Generated unique ID for this dataset
            user: Username/owner
            adata_path: Path to processed AnnData file
            tangram_adata_path: Path to tangram results (if applicable)
            geojson_path: Path or URL to GeoJSON
            genie_network_path: Path to GENIE3 network
            sponge_network_path: Path to SPONGE network

        Returns:
            VisiumDataset instance
        """
        if not adata_path:
            raise ValueError("adata_path is required to create a VisiumDataset")

        dataset = cls(
            dataset_id=dataset_id,
            alias=params.get_dataset_alias(),
            adata_path=adata_path,
            user=user,
            params=params,
            tangram_adata_path=tangram_adata_path,
            geojson_path=geojson_path,
            genie_network_path=genie_network_path,
            sponge_network_path=sponge_network_path,
            dataset_type=params.get_dataset_type(),
            use_tangram=params.is_tangram_enabled()
        )
        return dataset


# =============================================================================
# XeniumDataset - Xenium-Specific Implementation
# =============================================================================

class XeniumDataset(VisiumDataset):
    """
    Concrete implementation for Xenium spatial transcriptomics datasets.
    Extends VisiumDataset to handle Xenium-specific path selection:
    - When tangram is applied: uses adata_tg_scores.h5ad (tangram-projected data)
    - When tangram is not applied: uses xenium_cells_with_grid_scores.h5ad (grid-mapped cells)
    """

    def __init__(
        self,
        dataset_id: str,
        alias: str,
        adata_path: str,
        user: str = "anonymous",
        params: Optional[Params] = None,
        tangram_adata_path: Optional[str] = None,
        geojson_path: Optional[str] = None,
        genie_network_path: Optional[str] = None,
        sponge_network_path: Optional[str] = None,
        tf_graph_path: Optional[str] = None,
        created_at: Optional[datetime] = None,
        dataset_type: Optional[str] = None,
        use_tangram: bool = False,
        use_multiome: bool = False,
        xenium_grid_adata_path: Optional[str] = None,  # xenium_cells_with_grid_scores.h5ad
        **metadata
    ):
        """
        Initialize a XeniumDataset.

        Args:
            adata_path: Primary visualization adata path (set by use_tangram logic)
            xenium_grid_adata_path: Path to grid-mapped xenium cells (used when tangram=False)
            use_tangram: Whether tangram deconvolution was applied
            tf_graph_path: Optional path to TF graph JSON file
            **metadata: Additional Xenium-specific metadata

        Note:
            The primary adata_path is selected based on use_tangram:
            - If use_tangram=True: adata_path should point to adata_tg_scores.h5ad
            - If use_tangram=False: adata_path should point to xenium_cells_with_grid_scores.h5ad
              (or set via xenium_grid_adata_path and let from_params() decide)
        """
        super().__init__(
            dataset_id=dataset_id,
            alias=alias,
            adata_path=adata_path,
            user=user,
            params=params,
            tangram_adata_path=tangram_adata_path,
            geojson_path=geojson_path,
            genie_network_path=genie_network_path,
            sponge_network_path=sponge_network_path,
            tf_graph_path=tf_graph_path,
            created_at=created_at,
            dataset_type=dataset_type,
            use_tangram=use_tangram,
            use_multiome=use_multiome,
            **metadata
        )
        self.xenium_grid_adata_path = xenium_grid_adata_path

    def to_dict(self) -> Dict[str, Any]:
        """Serialize XeniumDataset to dictionary."""
        data = super().to_dict()
        data.update({
            "xenium_grid_adata_path": self.xenium_grid_adata_path,
        })
        return data

    @classmethod
    def from_params(
        cls,
        params: Params,
        dataset_id: str,
        user: str = "anonymous",
        adata_path: Optional[str] = None,
        tangram_adata_path: Optional[str] = None,
        geojson_path: Optional[str] = None,
        genie_network_path: Optional[str] = None,
        sponge_network_path: Optional[str] = None,
        xenium_grid_adata_path: Optional[str] = None,
        # Multiome-specific (accepted but ignored)
        adata_st_scores_path: Optional[str] = None,
        adata_tg_scores_path: Optional[str] = None,
    ) -> "XeniumDataset":
        """
        Create XeniumDataset from Params and output paths.

        Intelligently selects the primary adata_path based on tangram usage:
        - If tangram=True: uses tangram_adata_path (adata_tg_scores.h5ad)
        - If tangram=False: uses xenium_grid_adata_path (xenium_cells_with_grid_scores.h5ad)

        Args:
            params: Validated parameters
            dataset_id: Generated unique ID for this dataset
            user: Username/owner
            adata_path: If provided, use as primary path (overrides intelligent selection)
            tangram_adata_path: Path to tangram results (adata_tg_scores.h5ad)
            geojson_path: Path or URL to GeoJSON
            genie_network_path: Path to GENIE3 network
            sponge_network_path: Path to SPONGE network
            xenium_grid_adata_path: Path to grid-mapped xenium cells (xenium_cells_with_grid_scores.h5ad)

        Returns:
            XeniumDataset instance

        Raises:
            ValueError: If no valid adata path is provided
        """
        use_tangram_flag = params.is_tangram_enabled()

        # Determine primary adata path based on tangram usage
        if adata_path:
            # Explicit path provided, use it
            primary_adata_path = adata_path
        elif use_tangram_flag and tangram_adata_path:
            # Tangram was used, prefer tangram results
            primary_adata_path = tangram_adata_path
        elif xenium_grid_adata_path:
            # No tangram or tangram path not available, use grid adata
            primary_adata_path = xenium_grid_adata_path
        else:
            raise ValueError(
                "XeniumDataset requires adata_path, tangram_adata_path (if use_tangram=True), "
                "or xenium_grid_adata_path (if use_tangram=False)"
            )

        dataset = cls(
            dataset_id=dataset_id,
            alias=params.get_dataset_alias(),
            adata_path=primary_adata_path,
            user=user,
            params=params,
            tangram_adata_path=tangram_adata_path,
            geojson_path=geojson_path,
            genie_network_path=genie_network_path,
            sponge_network_path=sponge_network_path,
            dataset_type=params.get_dataset_type(),
            use_tangram=use_tangram_flag,
            xenium_grid_adata_path=xenium_grid_adata_path,
        )
        return dataset


# =============================================================================
# MultiomeDataset - Multiome (ATAC-seq) Implementation
# =============================================================================

class MultiomeDataset(Dataset):
    """
    Concrete implementation for Multiome (ATAC-seq) datasets.
    Extends Dataset with multiome-specific file paths from the pipeline outputs.

    Output files are written to: uploads/job_<id>/<outdir>/
    where <outdir> is a randomly generated directory name.
    """

    def __init__(
        self,
        dataset_id: str,
        alias: str,
        adata_path: str,
        user: str = "anonymous",
        params: Optional[Params] = None,
        tangram_adata_path: Optional[str] = None,
        geojson_path: Optional[str] = None,
        genie_network_path: Optional[str] = None,
        sponge_network_path: Optional[str] = None,
        created_at: Optional[datetime] = None,
        # Multiome-specific paths (always present)
        adata_st_scores_path: Optional[str] = None,
        adata_tg_scores_path: Optional[str] = None,
        adata_map_path: Optional[str] = None,
        adata_map_X_csv_path: Optional[str] = None,
        adata_map_var_csv_path: Optional[str] = None,
        calc_scores_log_path: Optional[str] = None,
        global_motif_analysis_path: Optional[str] = None,
        # Multiome-specific paths (conditional)
        motif_to_tf_csv_path: Optional[str] = None,
        spot_obj_chromvar_path: Optional[str] = None,
        spot_obj_footprints_path: Optional[str] = None,
        dissociated_obj_footprints_path: Optional[str] = None,
        chromvar_scores_csv_path: Optional[str] = None,
        diff_motif_activity_csv_paths: Optional[Dict[str, str]] = None,  # comparison -> path
        footprint_pdf_paths: Optional[Dict[str, str]] = None,  # motif_id -> path
        tf_graph_path: Optional[str] = None,
        grn_evaluation_name: Optional[str] = None,

        # Multiome flags
        use_chromvar: bool = False,
        use_differential_motif_activity: bool = False,
        use_footprinting: bool = False,
        use_moranI: bool = False,
        use_gearyC: bool = False,
        grn_evaluation_on_demand_available: bool = False,
        **metadata
    ):
        """
        Initialize a MultiomeDataset.

        Args:
            dataset_id: Unique identifier (e.g., "job_123_user")
            alias: User-friendly name
            adata_path: Path to main AnnData H5AD file (adata_st_scores.h5ad typically)
            user: User who created/owns this dataset
            params: Original Params object that created this dataset
            tangram_adata_path: Path to tangram-projected data (adata_tg_scores.h5ad)
            geojson_path: Optional path or URL to GeoJSON visualization
            genie_network_path: Optional path to GENIE3 network
            sponge_network_path: Optional path to SPONGE network
            created_at: Creation timestamp

            # Always-present multiome outputs
            adata_st_scores_path: AnnData for spatial transcriptomics with scores
            adata_tg_scores_path: AnnData for Tangram-projected data with scores
            adata_map_path: Raw Tangram mapping AnnData (cells × spots)
            adata_map_X_csv_path: Tangram mapping matrix as CSV
            adata_map_var_csv_path: Spatial spot metadata from Tangram
            calc_scores_log_path: Pipeline execution log
            global_motif_analysis_path: Complete Seurat multiome object (RDS)

            # Conditional multiome outputs
            motif_to_tf_csv_path: Motif ID to TF name mapping
            spot_obj_chromvar_path: Spot-level chromVAR results (RDS)
            spot_obj_footprints_path: Spot-level footprint profiles (RDS)
            dissociated_obj_footprints_path: Dissociated cell footprints (RDS)
            chromvar_scores_csv_path: chromVAR deviation scores (cells × motifs)
            diff_motif_activity_csv_paths: Dict mapping comparisons to CSV paths
            footprint_pdf_paths: Dict mapping motif IDs to footprint PDF paths
            tf_graph_path: Optional path to TF graph JSON file

            # Feature flags (which commands were executed in the pipeline)
            use_chromvar: Whether chromVAR analysis was performed
            use_differential_motif_activity: Whether differential motif activity was computed
            use_footprinting: Whether footprint analysis was performed
            use_moranI: Whether Moran's I spatial autocorrelation was computed
            use_gearyC: Whether Geary's C spatial autocorrelation was computed

            **metadata: Additional multiome-specific metadata
        """
        super().__init__(
            dataset_id=dataset_id,
            alias=alias,
            adata_path=adata_path,
            user=user,
            params=params,
            tangram_adata_path=tangram_adata_path,
            geojson_path=geojson_path,
            genie_network_path=genie_network_path,
            sponge_network_path=sponge_network_path,
            tf_graph_path=tf_graph_path,
            created_at=created_at,
            **metadata
        )

        # Always-present multiome outputs
        self.adata_st_scores_path = adata_st_scores_path
        self.adata_tg_scores_path = adata_tg_scores_path
        self.adata_map_path = adata_map_path
        self.adata_map_X_csv_path = adata_map_X_csv_path
        self.adata_map_var_csv_path = adata_map_var_csv_path
        self.calc_scores_log_path = calc_scores_log_path
        self.global_motif_analysis_path = global_motif_analysis_path

        # Conditional multiome outputs
        self.motif_to_tf_csv_path = motif_to_tf_csv_path
        self.spot_obj_chromvar_path = spot_obj_chromvar_path
        self.spot_obj_footprints_path = spot_obj_footprints_path
        self.dissociated_obj_footprints_path = dissociated_obj_footprints_path
        self.chromvar_scores_csv_path = chromvar_scores_csv_path
        self.diff_motif_activity_csv_paths = diff_motif_activity_csv_paths or {}
        self.footprint_pdf_paths = footprint_pdf_paths or {}

        # Feature flags indicating which pipeline steps were executed
        self.use_chromvar = use_chromvar
        self.use_differential_motif_activity = use_differential_motif_activity
        self.use_footprinting = use_footprinting
        self.use_moranI = use_moranI
        self.use_gearyC = use_gearyC
        self.grn_evaluation_on_demand_available = grn_evaluation_on_demand_available
        self.grn_evaluation_name = grn_evaluation_name

    def to_dict(self) -> Dict[str, Any]:
        """Serialize MultiomeDataset to dictionary."""
        data = super().to_dict()
        data.update({
            "adata_st_scores_path": self.adata_st_scores_path,
            "adata_tg_scores_path": self.adata_tg_scores_path,
            "adata_map_path": self.adata_map_path,
            "adata_map_X_csv_path": self.adata_map_X_csv_path,
            "adata_map_var_csv_path": self.adata_map_var_csv_path,
            "calc_scores_log_path": self.calc_scores_log_path,
            "global_motif_analysis_path": self.global_motif_analysis_path,
            "motif_to_tf_csv_path": self.motif_to_tf_csv_path,
            "spot_obj_chromvar_path": self.spot_obj_chromvar_path,
            "spot_obj_footprints_path": self.spot_obj_footprints_path,
            "dissociated_obj_footprints_path": self.dissociated_obj_footprints_path,
            "chromvar_scores_csv_path": self.chromvar_scores_csv_path,
            "diff_motif_activity_csv_paths": self.diff_motif_activity_csv_paths,
            "footprint_pdf_paths": self.footprint_pdf_paths,
            "tf_graph_path": self.tf_graph_path,
            "use_chromvar": self.use_chromvar,
            "use_differential_motif_activity": self.use_differential_motif_activity,
            "use_footprinting": self.use_footprinting,
            "use_moranI": self.use_moranI,
            "use_gearyC": self.use_gearyC,
            "grn_evaluation_on_demand_available": self.grn_evaluation_on_demand_available,
            "grn_evaluation_name": self.grn_evaluation_name,
        })
        return data

    def validate_paths(self) -> bool:
        """
        Verify that all referenced multiome files exist.
        Only checks paths that are set (not None).

        Returns:
            True if all existing paths are valid, False otherwise
        """
        # Always-required paths
        always_required = [
            ("adata_st_scores", self.adata_st_scores_path),
            ("adata_tg_scores", self.adata_tg_scores_path),
            ("adata_map", self.adata_map_path),
            ("adata_map.X.csv", self.adata_map_X_csv_path),
            ("adata_map.var.csv", self.adata_map_var_csv_path),
            ("global_motif_analysis", self.global_motif_analysis_path),
        ]

        all_valid = super().validate_paths()

        for name, path_str in always_required:
            if path_str:
                path = Path(path_str)
                if not path.exists():
                    print(f"✗ Dataset {self.id}: multiome file '{name}' not found at {path_str}")
                    all_valid = False

        # Optional paths - only warn if set but missing
        optional_paths = [
            ("motif_to_tf.csv", self.motif_to_tf_csv_path),
            ("spot_obj_chromvar", self.spot_obj_chromvar_path),
            ("spot_obj_footprints", self.spot_obj_footprints_path),
            ("dissociated_obj_footprints", self.dissociated_obj_footprints_path),
            ("chromvar_scores.csv", self.chromvar_scores_csv_path),
        ]

        for name, path_str in optional_paths:
            if path_str:
                path = Path(path_str)
                if not path.exists():
                    print(f"⚠ Dataset {self.id}: optional multiome file '{name}' not found at {path_str}")

        # Check dict-based paths
        for comparison, path_str in self.diff_motif_activity_csv_paths.items():
            if path_str:
                path = Path(path_str)
                if not path.exists():
                    print(f"⚠ Dataset {self.id}: diff motif activity CSV for '{comparison}' not found at {path_str}")

        for motif_id, path_str in self.footprint_pdf_paths.items():
            if path_str:
                path = Path(path_str)
                if not path.exists():
                    print(f"⚠ Dataset {self.id}: footprint PDF for motif '{motif_id}' not found at {path_str}")

        return all_valid

    @classmethod
    def from_params(
        cls,
        params: Params,
        dataset_id: str,
        user: str = "anonymous",
        adata_path: Optional[str] = None,
        tangram_adata_path: Optional[str] = None,
        adata_st_scores_path: Optional[str] = None,
        adata_tg_scores_path: Optional[str] = None,
        adata_map_path: Optional[str] = None,
        adata_map_X_csv_path: Optional[str] = None,
        adata_map_var_csv_path: Optional[str] = None,
        calc_scores_log_path: Optional[str] = None,
        global_motif_analysis_path: Optional[str] = None,
        motif_to_tf_csv_path: Optional[str] = None,
        spot_obj_chromvar_path: Optional[str] = None,
        spot_obj_footprints_path: Optional[str] = None,
        dissociated_obj_footprints_path: Optional[str] = None,
        chromvar_scores_csv_path: Optional[str] = None,
        diff_motif_activity_csv_paths: Optional[Dict[str, str]] = None,
        footprint_pdf_paths: Optional[Dict[str, str]] = None,
        geojson_path: Optional[str] = None,
        genie_network_path: Optional[str] = None,
        sponge_network_path: Optional[str] = None,
        xenium_grid_adata_path: Optional[str] = None,
    ) -> "MultiomeDataset":
        """
        Create MultiomeDataset from Params and multiome output paths.

        Args:
            params: Validated parameters
            dataset_id: Generated unique ID for this dataset
            user: Username/owner
            adata_path: Path to primary processed AnnData file (typically adata_st_scores.h5ad)
            tangram_adata_path: Alias for adata_tg_scores_path (ignored if adata_tg_scores_path is provided)
            adata_st_scores_path: Path to spatial transcriptomics AnnData with scores
            adata_tg_scores_path: Path to Tangram-projected AnnData with scores
            adata_map_path: Path to Tangram mapping AnnData
            adata_map_X_csv_path: Path to Tangram mapping matrix CSV
            adata_map_var_csv_path: Path to spatial spot metadata CSV
            calc_scores_log_path: Path to pipeline execution log
            global_motif_analysis_path: Path to Seurat multiome object (RDS)
            motif_to_tf_csv_path: Path to motif-to-TF mapping
            spot_obj_chromvar_path: Path to spot-level chromVAR object
            spot_obj_footprints_path: Path to spot-level footprints
            dissociated_obj_footprints_path: Path to dissociated cell footprints
            chromvar_scores_csv_path: Path to chromVAR scores CSV
            diff_motif_activity_csv_paths: Dict of comparison name -> CSV path
            footprint_pdf_paths: Dict of motif ID -> PDF path
            geojson_path: Path or URL to GeoJSON
            genie_network_path: Path to GENIE3 network
            sponge_network_path: Path to SPONGE network
            xenium_grid_adata_path: Xenium grid path (ignored for multiome datasets)

        Returns:
            MultiomeDataset instance
        """
        if not adata_st_scores_path and not adata_path:
            raise ValueError("Either adata_path or adata_st_scores_path is required to create a MultiomeDataset")

        # Use adata_path as fallback for adata_st_scores_path
        final_adata_st_scores = adata_st_scores_path or adata_path
        final_adata_path = adata_path or adata_st_scores_path

        # Use tangram_adata_path as fallback for adata_tg_scores_path (for compatibility)
        final_adata_tg_scores = adata_tg_scores_path or tangram_adata_path

        dataset = cls(
            dataset_id=dataset_id,
            alias=params.get_dataset_alias(),
            adata_path=final_adata_path,
            user=user,
            params=params,
            tangram_adata_path=final_adata_tg_scores,
            geojson_path=geojson_path,
            genie_network_path=genie_network_path,
            sponge_network_path=sponge_network_path,
            adata_st_scores_path=final_adata_st_scores,
            adata_tg_scores_path=final_adata_tg_scores,
            adata_map_path=adata_map_path,
            adata_map_X_csv_path=adata_map_X_csv_path,
            adata_map_var_csv_path=adata_map_var_csv_path,
            calc_scores_log_path=calc_scores_log_path,
            global_motif_analysis_path=global_motif_analysis_path,
            motif_to_tf_csv_path=motif_to_tf_csv_path,
            spot_obj_chromvar_path=spot_obj_chromvar_path,
            spot_obj_footprints_path=spot_obj_footprints_path,
            dissociated_obj_footprints_path=dissociated_obj_footprints_path,
            chromvar_scores_csv_path=chromvar_scores_csv_path,
            diff_motif_activity_csv_paths=diff_motif_activity_csv_paths,
            footprint_pdf_paths=footprint_pdf_paths,
            use_chromvar=bool(chromvar_scores_csv_path or spot_obj_chromvar_path),
            use_differential_motif_activity=bool(diff_motif_activity_csv_paths and len(diff_motif_activity_csv_paths) > 0),
            use_footprinting=bool(footprint_pdf_paths and len(footprint_pdf_paths) > 0),
            use_moranI=params.chromVar.moranI if params.chromVar else False,
            use_gearyC=params.chromVar.gearyC if params.chromVar else False,
        )
        return dataset

# =============================================================================
# DatasetFactory - Routing to Correct Dataset Type
# =============================================================================

class DatasetFactory:
    """
    Factory for creating Dataset objects of appropriate types.
    Routes based on dataset_type parameter.
    """

    DATASET_CLASSES = {
        "visium": VisiumDataset,
        "xenium": XeniumDataset,
        "multiome": MultiomeDataset,
    }

    @classmethod
    def create_dataset(
        cls,
        params: Params,
        dataset_id: str,
        user: str = "anonymous",
        adata_path: Optional[str] = None,
        tangram_adata_path: Optional[str] = None,
        geojson_path: Optional[str] = None,
        genie_network_path: Optional[str] = None,
        sponge_network_path: Optional[str] = None,
        # Xenium-specific paths
        xenium_grid_adata_path: Optional[str] = None,
        # Multiome-specific paths (always-present outputs)
        adata_st_scores_path: Optional[str] = None,
        adata_tg_scores_path: Optional[str] = None,
        adata_map_path: Optional[str] = None,
        adata_map_X_csv_path: Optional[str] = None,
        adata_map_var_csv_path: Optional[str] = None,
        calc_scores_log_path: Optional[str] = None,
        global_motif_analysis_path: Optional[str] = None,
        # Multiome-specific paths (conditional outputs)
        motif_to_tf_csv_path: Optional[str] = None,
        spot_obj_chromvar_path: Optional[str] = None,
        spot_obj_footprints_path: Optional[str] = None,
        dissociated_obj_footprints_path: Optional[str] = None,
        chromvar_scores_csv_path: Optional[str] = None,
        diff_motif_activity_csv_paths: Optional[Dict[str, str]] = None,
        footprint_pdf_paths: Optional[Dict[str, str]] = None,
    ) -> Dataset:
        """
        Create a Dataset of the appropriate type based on params.
        Routes to VisiumDataset, XeniumDataset, or MultiomeDataset based on dataset type.

        Args:
            params: Validated parameters containing dataset type info
            dataset_id: Unique identifier for the dataset
            user: Dataset owner/creator

            # Common paths
            adata_path: Path to processed AnnData file
            tangram_adata_path: Optional path to tangram output
            geojson_path: Optional geojson path/url
            genie_network_path: Optional GENIE3 network path
            sponge_network_path: Optional SPONGE network path

            # Xenium-specific
            xenium_grid_adata_path: Path to grid-mapped xenium cells (xenium_cells_with_grid_scores.h5ad)

            # Multiome: Always-present outputs
            adata_st_scores_path: Spatial transcriptomics AnnData with scores
            adata_tg_scores_path: Tangram-projected AnnData with scores
            adata_map_path: Tangram mapping AnnData (cells × spots)
            adata_map_X_csv_path: Tangram mapping matrix as CSV
            adata_map_var_csv_path: Spatial spot metadata CSV
            calc_scores_log_path: Pipeline execution log
            global_motif_analysis_path: Complete Seurat multiome object (RDS)

            # Multiome: Conditional outputs
            motif_to_tf_csv_path: Motif ID to TF name mapping
            spot_obj_chromvar_path: Spot-level chromVAR results (RDS)
            spot_obj_footprints_path: Spot-level footprints (RDS)
            dissociated_obj_footprints_path: Dissociated cell footprints (RDS)
            chromvar_scores_csv_path: chromVAR deviation scores (cells × motifs)
            diff_motif_activity_csv_paths: Dict mapping comparisons to CSV paths
            footprint_pdf_paths: Dict mapping motif IDs to footprint PDF paths

        Returns:
            Dataset instance (concrete subclass: VisiumDataset, XeniumDataset, or MultiomeDataset)

        Raises:
            ValueError: If dataset type is not supported
        """
        # Use registry type which considers multiome, but pass spatial type to params
        dataset_type = params.get_registry_dataset_type().lower()

        dataset_class = cls.DATASET_CLASSES.get(dataset_type, VisiumDataset)

        print(f"Creating {dataset_class.__name__} from params (type={dataset_type})")

        # Pass ALL parameters to from_params - each class handles only what it accepts
        dataset = dataset_class.from_params(
            params=params,
            dataset_id=dataset_id,
            user=user,
            adata_path=adata_path,
            tangram_adata_path=tangram_adata_path,
            geojson_path=geojson_path,
            genie_network_path=genie_network_path,
            sponge_network_path=sponge_network_path,
            # Xenium
            xenium_grid_adata_path=xenium_grid_adata_path,
            # Multiome always-present
            adata_st_scores_path=adata_st_scores_path,
            adata_tg_scores_path=adata_tg_scores_path,
            adata_map_path=adata_map_path,
            adata_map_X_csv_path=adata_map_X_csv_path,
            adata_map_var_csv_path=adata_map_var_csv_path,
            calc_scores_log_path=calc_scores_log_path,
            global_motif_analysis_path=global_motif_analysis_path,
            # Multiome conditional
            motif_to_tf_csv_path=motif_to_tf_csv_path,
            spot_obj_chromvar_path=spot_obj_chromvar_path,
            spot_obj_footprints_path=spot_obj_footprints_path,
            dissociated_obj_footprints_path=dissociated_obj_footprints_path,
            chromvar_scores_csv_path=chromvar_scores_csv_path,
            diff_motif_activity_csv_paths=diff_motif_activity_csv_paths,
            footprint_pdf_paths=footprint_pdf_paths,
        )

        return dataset

    @classmethod
    def from_registry_dict(cls, data: Dict[str, Any]) -> Dataset:
        """
        Deserialize from registry dictionary and route to appropriate type.

        Args:
            data: Dictionary from registry containing 'type' field

        Returns:
            Dataset instance of appropriate type
        """
        dataset_type = data.get("type", "visium").lower()
        dataset_class = cls.DATASET_CLASSES.get(dataset_type, VisiumDataset)
        return dataset_class.from_registry_dict(data)
