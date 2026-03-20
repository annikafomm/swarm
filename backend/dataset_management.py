import json
import os
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, List, Union, Any

from dataset_structure import Dataset, DatasetFactory, VisiumDataset, Params


class DatasetRegistry:
    """
    Manages available datasets across the application.
    Stores metadata about uploaded and built-in datasets.
    Holds Dataset objects in memory, persists as JSON to disk.
    """

    def __init__(self, registry_file: Path = None):
        if registry_file is None:
            # Use absolute path based on script location
            base_path = Path(__file__).parent
            registry_file = base_path / "uploads" / "dataset_registry.json"
        self.registry_file = registry_file
        self.registry_file.parent.mkdir(parents=True, exist_ok=True)
        print(f"Dataset registry path: {self.registry_file}")
        # self.datasets now stores Dataset objects organized by category
        self.datasets = self._load_registry()

    def _load_registry(self, clear_uploads: bool = True) -> Dict:
        """
        Load existing registry or create new one.
        Deserializes JSON into Dataset objects.

        Args:
            clear_uploads: If True (on startup), clear uploaded datasets.
                          If False (runtime), preserve existing uploads.

        Returns:
            Dict with "builtin" and "uploaded" keys containing Dataset objects
        """
        datasets = {"builtin": {}, "uploaded": {}}

        if self.registry_file.exists():
            try:
                with open(self.registry_file, 'r') as f:
                    data = json.load(f)

                    # Load builtin datasets
                    for dataset_id, dataset_dict in data.get("builtin", {}).items():
                        try:
                            # Mark as builtin to preserve category
                            dataset_dict["category"] = "builtin"
                            dataset = DatasetFactory.from_registry_dict(dataset_dict)
                            datasets["builtin"][dataset_id] = dataset
                        except Exception as e:
                            print(f"⚠ Failed to load builtin dataset {dataset_id}: {e}")

                    # Load uploaded datasets (if not clearing)
                    if not clear_uploads:
                        for dataset_id, dataset_dict in data.get("uploaded", {}).items():
                            try:
                                dataset_dict["category"] = "uploaded"
                                dataset = DatasetFactory.from_registry_dict(dataset_dict)
                                datasets["uploaded"][dataset_id] = dataset
                            except Exception as e:
                                print(f"⚠ Failed to load uploaded dataset {dataset_id}: {e}")

                return datasets

            except Exception as e:
                print(f"Error loading registry: {e}")
                return datasets

        return datasets

    def _save_registry(self):
        """
        Persist registry to disk.
        Converts Dataset objects to dictionaries for JSON serialization.
        """
        try:
            # Convert Dataset objects to dicts for JSON storage
            json_data = {
                "builtin": {
                    dataset_id: dataset.to_registry_dict()
                    for dataset_id, dataset in self.datasets.get("builtin", {}).items()
                },
                "uploaded": {
                    dataset_id: dataset.to_registry_dict()
                    for dataset_id, dataset in self.datasets.get("uploaded", {}).items()
                }
            }

            with open(self.registry_file, 'w') as f:
                json.dump(json_data, indent=2, fp=f, default=str)
            print(f"✓ Registry saved to {self.registry_file}")
        except Exception as e:
            print(f"✗ Error saving registry: {e}")
            raise

    def register_uploaded_dataset(
        self,
        dataset: Optional[Dataset] = None,
        dataset_id: Optional[str] = None,
        alias: Optional[str] = None,
        adata_path: Optional[str] = None,
        tangram_adata_path: Optional[str] = None,
        geojson_path: Optional[str] = None,
        genie_network_path: Optional[str] = None,
        sponge_network_path: Optional[str] = None,
        user: str = "anonymous",
        footprint_list: Optional[List[str]] = None,
        **metadata
    ) -> None:
        """
        Register a newly uploaded dataset.

        Can accept either:
        1. A Dataset object directly: register_uploaded_dataset(dataset=my_dataset)
        2. Legacy parameters: register_uploaded_dataset(dataset_id=..., alias=..., adata_path=..., etc.)

        Args:
            dataset: Dataset object to register (if provided, other args are ignored)
            dataset_id: Legacy: unique identifier
            alias: Legacy: user-friendly name
            adata_path: Legacy: path to AnnData H5AD file
            ... (other legacy parameters)
        """
        if "uploaded" not in self.datasets:
            self.datasets["uploaded"] = {}

        if dataset is not None:
            # Register Dataset object directly
            self.datasets["uploaded"][dataset.id] = dataset
            print(f"Registering Dataset object: {dataset.id}")
        else:
            # Legacy: create VisiumDataset from parameters
            if not all([dataset_id, alias, adata_path]):
                raise ValueError("Must provide either dataset object or (dataset_id, alias, adata_path)")

            dataset_obj = VisiumDataset(
                dataset_id=dataset_id,
                alias=alias,
                adata_path=adata_path,
                user=user,
                tangram_adata_path=tangram_adata_path,
                geojson_path=geojson_path,
                genie_network_path=genie_network_path,
                sponge_network_path=sponge_network_path,
                created_at=datetime.now(),
                footprint_list=footprint_list,
                **metadata
            )
            self.datasets["uploaded"][dataset_id] = dataset_obj
            print(f"Registered legacy dataset: {dataset_id}")

        print(f"Saving registry to {self.registry_file}")
        self._save_registry()
        print(f"✓ Registered dataset {dataset_id if dataset is None else dataset.id}")

    def register_builtin_dataset(
        self,
        dataset: Optional[Dataset] = None,
        dataset_id: Optional[str] = None,
        alias: Optional[str] = None,
        adata_path: Optional[str] = None,
        geojson_path: Optional[str] = None,
        genie_network_path: Optional[str] = None,
        sponge_network_path: Optional[str] = None,
        **metadata
    ) -> None:
        """
        Register a built-in dataset (from backend/data).

        Can accept either:
        1. A Dataset object directly: register_builtin_dataset(dataset=my_dataset)
        2. Legacy parameters: register_builtin_dataset(dataset_id=..., alias=..., adata_path=..., etc.)

        Args:
            dataset: Dataset object to register
            dataset_id: Legacy: unique identifier
            alias: Legacy: user-friendly name
            adata_path: Legacy: path to AnnData H5AD file
            ... (other legacy parameters)
        """
        if "builtin" not in self.datasets:
            self.datasets["builtin"] = {}

        if dataset is not None:
            # Register Dataset object directly
            self.datasets["builtin"][dataset.id] = dataset
            print(f"Registering builtin Dataset: {dataset.id}")
        else:
            # Legacy: create VisiumDataset from parameters
            if not all([dataset_id, alias, adata_path]):
                raise ValueError("Must provide either dataset object or (dataset_id, alias, adata_path)")

            dataset_obj = VisiumDataset(
                dataset_id=dataset_id,
                alias=alias,
                adata_path=adata_path,
                user="builtin",
                geojson_path=geojson_path,
                genie_network_path=genie_network_path,
                sponge_network_path=sponge_network_path,
                created_at=datetime.now(),
                **metadata
            )
            self.datasets["builtin"][dataset_id] = dataset_obj
            print(f"Registered legacy builtin dataset: {dataset_id}")

        self._save_registry()

    def get_all_datasets(self, as_dict: bool = True) -> Dict:
        """
        Return all available datasets grouped by type.

        Args:
            as_dict: If True, return Dataset objects serialized as dicts.
                    If False, return Dataset objects directly.

        Returns:
            Dict: {"builtin": {...}, "uploaded": {...}} with Dataset entries
        """
        if not as_dict:
            # Return Dataset objects directly (for in-app use)
            return self.datasets

        # Return serialized dicts (for JSON API responses)
        return {
            "builtin": {
                dataset_id: dataset.to_dict()
                for dataset_id, dataset in self.datasets.get("builtin", {}).items()
            },
            "uploaded": {
                dataset_id: dataset.to_dict()
                for dataset_id, dataset in self.datasets.get("uploaded", {}).items()
            }
        }

    def get_dataset_by_id(self, dataset_id: str, as_dict: bool = True) -> Optional[Union[Dataset, Dict]]:
        """
        Get a single dataset by ID.

        Args:
            dataset_id: Dataset identifier
            as_dict: If True, return serialized dict. If False, return Dataset object.

        Returns:
            Dataset object (if as_dict=False) or dict (if as_dict=True), or None if not found
        """
        for category in ["builtin", "uploaded"]:
            if category in self.datasets and dataset_id in self.datasets[category]:
                dataset = self.datasets[category][dataset_id]
                return dataset.to_dict() if as_dict else dataset
        return None

    def get_user_datasets(self, user: str, as_dict: bool = True) -> Dict[str, Union[Dataset, Dict]]:
        """
        Get all datasets for a specific user.

        Args:
            user: Username to filter by
            as_dict: If True, return serialized dicts. If False, return Dataset objects.

        Returns:
            Dict of {dataset_id: Dataset or dict} for datasets owned by user
        """
        user_datasets = {}
        if "uploaded" in self.datasets:
            for dataset_id, dataset in self.datasets["uploaded"].items():
                if dataset.user == user:
                    user_datasets[dataset_id] = dataset.to_dict() if as_dict else dataset
        return user_datasets

    def remove_dataset(self, dataset_id: str) -> bool:
        """Remove a dataset from registry"""
        for category in ["builtin", "uploaded"]:
            if category in self.datasets and dataset_id in self.datasets[category]:
                del self.datasets[category][dataset_id]
                self._save_registry()
                return True
        return False

    def update_dataset_paths(self, dataset_id: str, **new_paths) -> bool:
        """
        Update paths for an existing dataset.

        Args:
            dataset_id: Dataset to update
            **new_paths: Keyword arguments for Dataset attributes to update
                        (e.g., adata_path=..., geojson_path=...)

        Returns:
            True if update successful, False if dataset not found
        """
        dataset = self.get_dataset_by_id(dataset_id, as_dict=False)
        if dataset:
            for key, value in new_paths.items():
                if hasattr(dataset, key):
                    setattr(dataset, key, value)
                else:
                    print(f"⚠ Dataset attribute '{key}' does not exist")
            self._save_registry()
            return True
        return False

    def find_config_files(self, uploads_dir: Path) -> List[Path]:
        """
        Scan uploads directory for unregistered configuration files.
        Returns list of config files for datasets not yet in registry.

        Args:
            uploads_dir: Path to uploads directory

        Returns:
            List of Paths to unregistered config files
        """
        if not uploads_dir.exists():
            print(f"⚠ Uploads directory does not exist: {uploads_dir}")
            return []

        unregistered = []
        uploads_dir = Path(uploads_dir)

        print(f"[DEBUG] Scanning for config files in: {uploads_dir}")
        print(f"[DEBUG] Currently registered uploaded datasets: {list(self.datasets.get('uploaded', {}).keys())}")

        # Find all config files matching pattern job_*_config.json
        all_config_files = list(uploads_dir.rglob("job_*_config.json"))
        print(f"[DEBUG] Found {len(all_config_files)} config files total")

        for config_file in all_config_files:
            try:
                # Extract dataset_id from filename and directory
                # File: job_TIMESTAMP_UUID/job_TIMESTAMP_config.json
                parent_dir = config_file.parent
                dataset_id = parent_dir.name

                print(f"[DEBUG] Checking config: {dataset_id} - Registered: {dataset_id in self.datasets.get('uploaded', {})}")

                # Skip if already registered
                if dataset_id in self.datasets.get("uploaded", {}):
                    print(f"[DEBUG] Skipping {dataset_id} - already registered")
                    continue

                unregistered.append(config_file)
                print(f"[DEBUG] Added unregistered: {dataset_id}")
            except Exception as e:
                print(f"⚠ Error processing config file {config_file}: {e}")

        print(f"[DEBUG] Total unregistered config files found: {len(unregistered)}")
        return unregistered

    def get_unregistered_datasets(self, uploads_dir: Path) -> List[Dict[str, Any]]:
        """
        Get summary info about unregistered datasets from config files.
        Suitable for displaying in UI popup.

        Args:
            uploads_dir: Path to uploads directory

        Returns:
            List of dicts with dataset info for UI display
        """
        config_files = self.find_config_files(uploads_dir)
        unregistered = []

        for config_file in config_files:
            try:
                with open(config_file, 'r') as f:
                    config = json.load(f)

                parent_dir = config_file.parent
                dataset_id = parent_dir.name

                # Extract key info
                dataset_type = config.get("dataset", "Unknown")
                adata_path = config.get("output_files", {}).get("adata_path")
                created_at = datetime.fromtimestamp(
                    int(dataset_id.split("_")[1]) / 1000, # timestamp in milliseconds
                ) if "_" in dataset_id else datetime.now()

                # Check if adata file exists
                adata_exists = Path(adata_path).exists() if adata_path else False

                unregistered.append({
                    "dataset_id": dataset_id,
                    "alias": f"{dataset_type} Dataset - {dataset_id}",
                    "dataset_type": dataset_type,
                    "created_at": created_at.isoformat(),
                    "config_path": str(config_file),
                    "adata_path": adata_path,
                    "adata_exists": adata_exists,
                    "status": "ready" if adata_exists else "missing_files",
                })
            except Exception as e:
                print(f"⚠ Error reading unregistered dataset from {config_file}: {e}")

        # Sort by creation date, newest first
        unregistered.sort(key=lambda x: x["created_at"], reverse=True)
        return unregistered

    def register_dataset_from_config(
        self,
        config_file: Path,
        user: str = "anonymous",
    ) -> Dataset:
        """
        Register a dataset by loading its configuration file.
        Creates a Dataset object and adds it to the registry.

        Args:
            config_file: Path to config JSON file (e.g., job_123_config.json)
            user: Username to assign as owner

        Returns:
            Registered Dataset object

        Raises:
            ValueError: If config file is invalid or adata file missing
        """
        config_file = Path(config_file)
        if not config_file.exists():
            raise ValueError(f"Config file not found: {config_file}")

        with open(config_file, 'r') as f:
            config = json.load(f)

        parent_dir = config_file.parent
        dataset_id = parent_dir.name

        # Check if already registered
        if dataset_id in self.datasets.get("uploaded", {}):
            raise ValueError(f"Dataset {dataset_id} already registered")

        # Extract paths and metadata from config
        adata_path = config.get("output_files", {}).get("adata_path")
        geojson_path = config.get("output_files", {}).get("geojson_path")
        genie_network_path = config.get("output_files", {}).get("genie_network_path")
        sponge_network_path = config.get("output_files", {}).get("sponge_network_path")

        # Validate adata file exists
        if not adata_path or not Path(adata_path).exists():
            raise ValueError(f"AnnData file not found: {adata_path}")

        # Extract processing options
        dataset_type = config.get("dataset", "Visium")
        use_tangram = config.get("tangram", {}).get("use", False) if config.get("tangram") else False
        use_multiome = config.get("multiome", {}).get("use", False) if config.get("multiome") else False

        # Create Dataset object
        dataset = VisiumDataset(
            dataset_id=dataset_id,
            alias=f"{dataset_type} Dataset",
            adata_path=adata_path,
            user=user,
            tangram_adata_path=None,  # Can be updated if available
            geojson_path=geojson_path,
            genie_network_path=genie_network_path,
            sponge_network_path=sponge_network_path,
            created_at=datetime.now(),
            dataset_type=dataset_type,
            use_tangram=use_tangram,
            use_multiome=use_multiome,
        )

        # Register it
        self.register_uploaded_dataset(dataset=dataset)
        print(f"✓ Registered dataset from config: {dataset_id}")
        return dataset

    def delete_unregistered_dataset(self, dataset_id: str, uploads_dir: Path, delete_files: bool = True) -> bool:
        """
        Delete an unregistered dataset's directory and optionally its files.

        Args:
            dataset_id: Dataset identifier (directory name)
            uploads_dir: Path to uploads directory
            delete_files: If True, delete the entire job directory. If False, just remove references.

        Returns:
            True if deletion successful, False otherwise
        """
        try:
            job_dir = Path(uploads_dir) / dataset_id

            if not job_dir.exists():
                print(f"⚠ Dataset directory not found: {job_dir}")
                return False

            if delete_files:
                import shutil
                shutil.rmtree(job_dir)
                print(f"✓ Deleted dataset directory: {job_dir}")
            else:
                print(f"⚠ Dataset directory preserved: {job_dir}")

            return True
        except Exception as e:
            print(f"✗ Error deleting dataset {dataset_id}: {e}")
            return False

    def rescan_uploads_folder(self, uploads_dir: Path) -> Dict[str, str]:
        """
        TESTING ONLY: Scan the uploads folder and re-register datasets found there.
        Looks for *_st_scores.h5ad files and attempts to reconstruct metadata.
        Creates Dataset objects from found files.

        Returns a dict with results: {dataset_id: status_message}
        """
        if not uploads_dir.exists():
            return {"error": f"Uploads directory not found: {uploads_dir}"}

        results = {}
        uploads_dir = Path(uploads_dir)

        # Find all h5ad files
        for h5ad_file in uploads_dir.rglob("*_st_scores.h5ad"):
            try:
                # Extract dataset_id from parent directory name
                # Format: job_TIMESTAMP_USER/...filename
                job_dir = h5ad_file.parent.parent
                dataset_id = job_dir.name  # e.g., "job_123456_user"

                if dataset_id in self.datasets.get("uploaded", {}):
                    results[dataset_id] = "already registered"
                    continue

                # Try to load corresponding config file
                config_file = job_dir / f"{dataset_id}_config.json"
                metadata = {}
                geojson_path = None
                genie_network_path = None
                sponge_network_path = None
                dataset_type = "visium"
                use_tangram = False
                use_multiome = False

                if config_file.exists():
                    try:
                        with open(config_file, 'r') as f:
                            config = json.load(f)
                            dataset_type = config.get("dataset", "visium")
                            use_tangram = config.get("tangram", {}).get("use", False) if config.get("tangram") else False
                            use_multiome = config.get("multiome", {}).get("use", False) if config.get("multiome") else False
                            if "output_files" in config:
                                genie_network_path = config["output_files"].get("genie_network_path")
                                sponge_network_path = config["output_files"].get("sponge_network_path")
                    except Exception as e:
                        print(f"  ⚠ Could not parse config for {dataset_id}: {e}")

                # Look for geojson in multiple locations
                search_dir = h5ad_file.parent
                max_depth = 5
                for _ in range(max_depth):
                    candidates = list(search_dir.glob("hexagons.geojson"))
                    if candidates:
                        geojson_path = f"/api/geojson/{dataset_id}"
                        break
                    if search_dir <= job_dir:
                        break
                    search_dir = search_dir.parent

                # Extract user from directory name
                user_part = dataset_id.split("_", 2)[-1] if "_" in dataset_id else "unknown"

                # Create Dataset object
                dataset = VisiumDataset(
                    dataset_id=dataset_id,
                    alias=f"[RESCANNED] {dataset_id}",
                    adata_path=str(h5ad_file),
                    user=user_part,
                    geojson_path=geojson_path,
                    genie_network_path=genie_network_path,
                    sponge_network_path=sponge_network_path,
                    created_at=datetime.now(),
                    dataset_type=dataset_type,
                    use_tangram=use_tangram,
                    use_multiome=use_multiome,
                )

                # Register the Dataset object
                self.register_uploaded_dataset(dataset=dataset)
                results[dataset_id] = "registered"

            except Exception as e:
                results[str(h5ad_file)] = f"error: {e}"
                import traceback
                traceback.print_exc()

        return results
