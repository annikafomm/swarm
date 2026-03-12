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
        """Extract dataset type string."""
        return self.dataset.value if isinstance(self.dataset, Enum) else str(self.dataset)

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
            "created_at": self.created_at.isoformat(),
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
        Routes to appropriate subclass based on 'type' field.

        Args:
            data: Dictionary from registry (should have 'type' field)

        Returns:
            Appropriate Dataset subclass instance
        """
        data_copy = data.copy()
        dataset_type = data_copy.pop("type", "visium")

        # Get dataset class based on type
        # Currently all types use VisiumDataset, but structure supports inheritance
        from dataset_structure import DatasetFactory
        dataset_class = VisiumDataset

        # Handle 'id' vs 'dataset_id' field name
        if "id" in data_copy and "dataset_id" not in data_copy:
            data_copy["dataset_id"] = data_copy.pop("id")

        return dataset_class(**data_copy)

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
        created_at: Optional[datetime] = None,
        dataset_type: Optional[str] = None,  # e.g., "Visium", "Xenium"
        use_tangram: bool = False,
        use_multiome: bool = False,
        **metadata
    ):
        """
        Initialize a VisiumDataset.

        Args:
            dataset_type: Original dataset type from params (e.g., "Visium", "XeniumRIF1")
            use_tangram: Whether tangram deconvolution was applied
            use_multiome: Whether multiome data was processed
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
            created_at=created_at,
            **metadata
        )
        self.dataset_type = dataset_type  # e.g., "Visium"
        self.use_tangram = use_tangram
        self.use_multiome = use_multiome

    def to_dict(self) -> Dict[str, Any]:
        """Serialize VisiumDataset to dictionary."""
        data = super().to_dict()
        data.update({
            "dataset_type": self.dataset_type,
            "use_tangram": self.use_tangram,
            "use_multiome": self.use_multiome,
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
            use_tangram=params.is_tangram_enabled(),
            use_multiome=params.is_multiome_enabled(),
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
        "xenium": VisiumDataset,  # Placeholder: will create XeniumDataset subclass later
        "multiome": VisiumDataset,  # Placeholder: will create MultiomeDataset subclass later
        "singlecell": VisiumDataset,  # Placeholder: will create SingleCellDataset subclass later
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
    ) -> Dataset:
        """
        Create a Dataset of the appropriate type based on params.

        Args:
            params: Validated parameters containing dataset type info
            dataset_id: Unique identifier for the dataset
            user: Dataset owner/creator
            adata_path: Path to processed AnnData file
            tangram_adata_path: Optional path to tangram output
            geojson_path: Optional geojson path/url
            genie_network_path: Optional GENIE3 network path
            sponge_network_path: Optional SPONGE network path

        Returns:
            Dataset instance (concrete subclass)

        Raises:
            ValueError: If dataset type is not supported
        """
        dataset_type = params.get_dataset_type().lower()

        dataset_class = cls.DATASET_CLASSES.get(dataset_type, VisiumDataset)

        print(f"Creating {dataset_class.__name__} from params (type={dataset_type})")

        dataset = dataset_class.from_params(
            params=params,
            dataset_id=dataset_id,
            user=user,
            adata_path=adata_path,
            tangram_adata_path=tangram_adata_path,
            geojson_path=geojson_path,
            genie_network_path=genie_network_path,
            sponge_network_path=sponge_network_path,
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
