import json
import os
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, List, Union, Any

from dataset_structure import Dataset, DatasetFactory, VisiumDataset, XeniumDataset, MultiomeDataset, Params


class DatasetRegistry:
    """
    Manages available datasets across the application.
    Stores metadata about uploaded and built-in datasets.
    Holds Dataset objects in memory, persists as JSON to disk.
    """

    def __init__(self, registry_file: Path = None):
        # Disk persistence (dataset_registry.json) is disabled — see _load_registry
        # and _save_registry below. Kept commented out rather than deleted since the
        # registry_file plumbing would need to come back if that ever changes.
        # if registry_file is None:
        #     # Use absolute path based on script location
        #     base_path = Path(__file__).parent
        #     registry_file = base_path / "uploads" / "dataset_registry.json"
        # self.registry_file = registry_file
        # self.registry_file.parent.mkdir(parents=True, exist_ok=True)
        # print(f"Dataset registry path: {self.registry_file}")
        # self.datasets now stores Dataset objects organized by category
        self.datasets = self._load_registry()

    def _load_registry(self, clear_uploads: bool = True) -> Dict:
        """
        Load existing registry or create new one.
        Deserializes JSON into Dataset objects.

        Disk persistence is disabled (see below) — this always returns empty
        buckets. Builtins are rebuilt fresh from DEFAULT_DATASETS + their
        config/files by main.py's startup lifespan hook; uploaded datasets are
        recovered from job_*_config.json via find_config_files/register_dataset_from_config
        instead, which also verifies the referenced output files still exist on disk
        (dataset_registry.json never did — it was a write-only snapshot nothing read back,
        since this method was always called with clear_uploads=True).

        Args:
            clear_uploads: If True (on startup), clear uploaded datasets.
                          If False (runtime), preserve existing uploads.

        Returns:
            Dict with "builtin" and "uploaded" keys containing Dataset objects
            (both always start empty; "builtin" is populated later by the startup hook)
        """
        datasets = {"builtin": {}, "uploaded": {}}

        # if self.registry_file.exists():
        #     try:
        #         with open(self.registry_file, 'r') as f:
        #             data = json.load(f)
        #
        #             # Load uploaded datasets (if not clearing)
        #             if not clear_uploads:
        #                 for dataset_id, dataset_dict in data.get("uploaded", {}).items():
        #                     try:
        #                         dataset_dict["category"] = "uploaded"
        #                         dataset = DatasetFactory.from_registry_dict(dataset_dict)
        #                         datasets["uploaded"][dataset_id] = dataset
        #                     except Exception as e:
        #                         print(f"⚠ Failed to load uploaded dataset {dataset_id}: {e}")
        #
        #         return datasets
        #
        #     except Exception as e:
        #         print(f"Error loading registry: {e}")
        #         return datasets

        return datasets

    def _save_registry(self):
        """
        Disk persistence is disabled — this is now a no-op.
        Previously converted Dataset objects to dicts and wrote them to
        dataset_registry.json, but nothing ever read that file back (_load_registry
        is always called with clear_uploads=True), so it was pure write overhead.
        Kept commented out rather than deleted in case disk persistence is revisited.
        """
        pass
        # try:
        #     # Convert Dataset objects to dicts for JSON storage
        #     json_data = {
        #         "uploaded": {
        #             dataset_id: dataset.to_registry_dict()
        #             for dataset_id, dataset in self.datasets.get("uploaded", {}).items()
        #         }
        #     }
        #
        #     with open(self.registry_file, 'w') as f:
        #         json.dump(json_data, indent=2, fp=f, default=str)
        #     print(f"✓ Registry saved to {self.registry_file}")
        # except Exception as e:
        #     print(f"✗ Error saving registry: {e}")
        #     raise

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
            print()

        # print(f"Saving registry to {self.registry_file}")
        # self._save_registry()
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

        # self._save_registry()

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
                # self._save_registry()
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
            # self._save_registry()
            return True
        return False

    def find_config_files(self, uploads_dir: Path, user: Optional[str] = None) -> List[Path]:
        """
        Scan uploads directory for unregistered configuration files.
        Returns list of config files for datasets not yet in registry.

        Args:
            uploads_dir: Path to uploads directory
            user: Optional username to filter datasets for

        Returns:
            List of Paths to unregistered config files
        """
        if not uploads_dir.exists():
            print(f"⚠ Uploads directory does not exist: {uploads_dir}")
            return []

        unregistered = []
        uploads_dir = Path(uploads_dir)

        print(f"[DEBUG] Scanning for config files in: {uploads_dir} (filter user: {user})")
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

                # Filter by user if specified
                if user is not None:
                    # Check if directory name ends with _{user} or matches user
                    if not (dataset_id.endswith(f"_{user}") or dataset_id == user):
                        # Also check config file content if user field is present
                        try:
                            with open(config_file, 'r') as f:
                                cfg = json.load(f)
                                cfg_user = cfg.get("user") or cfg.get("email")
                                if cfg_user != user:
                                    continue
                        except Exception:
                            continue

                print(f"[DEBUG] Checking config: {dataset_id} - Registered: {dataset_id in self.datasets.get('uploaded', {})}")

                # Allow re-registration even if already registered
                # This enables recovery and re-registration on new sessions
                unregistered.append(config_file)
                print(f"[DEBUG] Added dataset for registration: {dataset_id}")
            except Exception as e:
                print(f"⚠ Error processing config file {config_file}: {e}")

        print(f"[DEBUG] Total unregistered config files found: {len(unregistered)}")
        return unregistered

    def get_unregistered_datasets(self, uploads_dir: Path, user: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        Get summary info about unregistered datasets from config files.
        Suitable for displaying in UI popup.
        Intelligently checks for the correct adata file based on dataset type and processing flags.
        Validates multiome and xenium-specific required files.

        Args:
            uploads_dir: Path to uploads directory
            user: Optional username to filter datasets for

        Returns:
            List of dicts with dataset info for UI display
        """
        config_files = self.find_config_files(uploads_dir, user=user)
        unregistered = []

        for config_file in config_files:
            try:
                with open(config_file, 'r') as f:
                    config = json.load(f)

                parent_dir = config_file.parent
                dataset_id = parent_dir.name

                # Extract key info
                dataset_type = config.get("dataset", "Unknown")
                if dataset_type == "Visium" and config.get("multiome").get("use") == True:
                    dataset_type = "Multiome"
                output_files = config.get("output_files", {})


                # Determine which adata file to check based on dataset type and processing
                adata_path = self._determine_adata_path(dataset_type, output_files, config)

                created_at = datetime.fromtimestamp(
                    int(dataset_id.split("_")[1]) / 1000, # timestamp in milliseconds
                ) if "_" in dataset_id else datetime.now()

                # Check if adata file exists
                adata_exists = Path(adata_path).exists() if adata_path else False
                print(f"[DEBUG] Dataset {dataset_id} - Type: {dataset_type} - Checking adata at: {adata_path} - Exists: {adata_exists}")

                # Validate type-specific required files
                missing_files = []
                is_complete = adata_exists
                dataset_type_lower = dataset_type.lower()

                if "multiome" in dataset_type_lower:
                    # Check multiome-specific paths
                    st_scores_path = output_files.get("adata_st_scores_path")
                    tg_scores_path = output_files.get("adata_tg_scores_path")
                    if not st_scores_path and not tg_scores_path:
                        missing_files.append("adata_st_scores.h5ad or adata_tg_scores.h5ad")
                        is_complete = False
                elif "xenium" in dataset_type_lower:
                    # Check xenium-specific paths
                    xenium_grid_path = output_files.get("xenium_grid_adata_path")
                    xenium_grid_exists = Path(xenium_grid_path).exists() if xenium_grid_path else False
                    if not xenium_grid_exists and not adata_exists:
                        missing_files.append("xenium_cells_with_grid_scores.h5ad or adata_path")
                        is_complete = False
                # For visium, just check main adata
                elif "visium" in dataset_type_lower:
                    missing_files.append("adata.h5ad")
                print(f"[DEBUG] Dataset {dataset_id} - Missing files: {missing_files} - Is complete: {is_complete}")

                status = "ready" if is_complete else "missing_files"
                missing_files_str = ", ".join(missing_files) if missing_files else None

                unregistered.append({
                    "dataset_id": dataset_id,
                    "alias": f"{dataset_type} Dataset - {dataset_id}",
                    "dataset_type": dataset_type,
                    "created_at": created_at.isoformat(),
                    "config_path": str(config_file),
                    "adata_path": adata_path,
                    "adata_exists": adata_exists,
                    "status": status,
                    "missing_files": missing_files_str,
                    "is_complete": is_complete,
                })
            except Exception as e:
                print(f"⚠ Error reading unregistered dataset from {config_file}: {e}")

        # Sort by creation date, newest first
        unregistered.sort(key=lambda x: x["created_at"], reverse=True)
        return unregistered

    def _determine_adata_path(self, dataset_type: str, output_files: Dict[str, Any], config: Dict[str, Any]) -> Optional[str]:
        """
        Determine which adata file to check based on dataset type and processing flags.
        (For registering datasets from config)
        Args:
            dataset_type: Type of dataset (e.g., "Visium", "Xenium", "Multiome")
            output_files: Dict of output file paths from config
            config: Full config dict (for processing flags like use_tangram)

        Returns:
            Path to the primary adata file to validate
        """
        # Handle Xenium datasets - they use different adata files based on tangram
        if "xenium" in dataset_type.lower():
            use_tangram = config.get("tangram", {}).get("use", False) if config.get("tangram") else False

            if use_tangram:
                # Tangram was applied, check for adata_tg_scores.h5ad
                return output_files.get("adata_tg_scores_path")
            else:
                # No tangram, check for xenium_cells_with_grid_scores.h5ad
                # This might be stored as adata_path or xenium_grid_adata_path in config
                return output_files.get("xenium_grid_adata_path") or output_files.get("adata_path")

        # Handle Multiome datasets - use adata_st_scores_path if available
        if "multiome" in dataset_type.lower():
            print(f"Multiome dataset detected [DEBUG], checking for adata_st_scores_path and adata_tg_scores_path")
            return output_files.get("adata_tg_scores_path")

        # Handle Visium datasets - prefer adata_st_scores_path (original ST with scores) if available
        if "visium" in dataset_type.lower():
            use_tangram = config.get("tangram", {}).get("use", False) if config.get("tangram") else False

            if use_tangram:
                # Tangram was applied, prefer adata_tg_scores (tangram-projected with scores)
                if output_files.get("adata_tg_scores_path"):
                    return output_files.get("adata_tg_scores_path")
            else:
                # No tangram, prefer adata_st_scores (original ST with scores)
                if output_files.get("adata_st_scores_path"):
                    return output_files.get("adata_st_scores_path")

            # Fallback to original adata if scores not available
            return output_files.get("adata_path")

        # Default for other types
        return output_files.get("adata_path")

    def register_dataset_from_config(
        self,
        config_file: Path,
        user: str = "anonymous",
    ) -> Dataset:
        """
        Register a dataset by loading its configuration file.
        Creates a Dataset object of the appropriate type and adds it to the registry.

        Args:
            config_file: Path to config JSON file (e.g., job_123_config.json)
            user: Username to assign as owner

        Returns:
            Registered Dataset object (concrete type varies: VisiumDataset, XeniumDataset, etc.)

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

        # Allow re-registration of already registered datasets
        # This supports recovery flows and re-registration on new sessions
        if dataset_id in self.datasets.get("uploaded", {}):
            print(f"[DEBUG] Dataset {dataset_id} already registered, re-registering...")

        # Extract paths and metadata from config
        output_files = config.get("output_files", {})
        dataset_type = config.get("dataset", "Visium")

        if dataset_type == "Visium" and config.get("multiome", {}).get("use") == True:
            dataset_type = "Multiome"

        # Determine which adata file to use based on dataset type
        adata_path = self._determine_adata_path(dataset_type, output_files, config)

        # For xenium datasets, also extract grid path for XeniumDataset
        xenium_grid_adata_path = output_files.get("xenium_grid_adata_path") or output_files.get("adata_path")

        geojson_path = output_files.get("geojson_path")
        genie_network_path = output_files.get("genie_network_path")
        sponge_network_path = output_files.get("sponge_network_path")

        # Validate adata file exists
        if not adata_path or not Path(adata_path).exists():
            raise ValueError(f"AnnData file not found: {adata_path}")

        # Extract processing options
        use_tangram = config.get("tangram", {}).get("use", False) if config.get("tangram") else False
        use_multiome = config.get("multiome", {}).get("use", False) if config.get("multiome") else False

        # Use DatasetFactory to create the appropriate dataset type
        # First, create a minimal Params object for the factory
        # (The factory expects Params but we'll pass what we can from config)

        # Route based on dataset type
        dataset_type_lower = dataset_type.lower()

        if "xenium" in dataset_type_lower:
            # Create XeniumDataset with both adata paths
            dataset = XeniumDataset(
                dataset_id=dataset_id,
                alias=f"{dataset_type} Dataset",
                adata_path=adata_path,
                user=user,
                tangram_adata_path=output_files.get("adata_tg_scores_path") if use_tangram else None,
                geojson_path=geojson_path,
                genie_network_path=genie_network_path,
                sponge_network_path=sponge_network_path,
                tf_graph_path=output_files.get("tf_graph_path"),
                created_at=datetime.now(),
                dataset_type=dataset_type,
                use_tangram=use_tangram,
                use_multiome=use_multiome,
                xenium_grid_adata_path=xenium_grid_adata_path,
            )
        elif "multiome" in dataset_type_lower:
            # Create MultiomeDataset with all multiome paths
            dataset = MultiomeDataset(
                dataset_id=dataset_id,
                alias=f"{dataset_type} Dataset",
                adata_path=adata_path,
                user=user,
                tangram_adata_path=output_files.get("adata_tg_scores_path"),
                geojson_path=geojson_path,
                genie_network_path=genie_network_path,
                sponge_network_path=sponge_network_path,
                created_at=datetime.now(),
                adata_st_scores_path=output_files.get("adata_st_scores_path"),
                adata_tg_scores_path=output_files.get("adata_tg_scores_path"),
                adata_map_path=output_files.get("adata_map_path"),
                adata_map_X_csv_path=output_files.get("adata_map_X_csv_path"),
                adata_map_var_csv_path=output_files.get("adata_map_var_csv_path"),
                calc_scores_log_path=output_files.get("calc_scores_log_path"),
                global_motif_analysis_path=output_files.get("global_motif_analysis_path"),
                motif_to_tf_csv_path=output_files.get("motif_to_tf_csv_path"),
                spot_obj_chromvar_path=output_files.get("spot_obj_chromvar_path"),
                spot_obj_footprints_path=output_files.get("spot_obj_footprints_path"),
                dissociated_obj_footprints_path=output_files.get("dissociated_obj_footprints_path"),
                chromvar_scores_csv_path=output_files.get("chromvar_scores_csv_path"),
                diff_motif_activity_csv_paths=output_files.get("diff_motif_activity_csv_paths"),
                footprint_pdf_paths=output_files.get("footprint_pdf_paths"),
                tf_graph_path=output_files.get("tf_graph_path"),
                grn_evaluation_name = output_files.get("grn_evaluation_name")
            )
        else:
            # Default to VisiumDataset
            # For Visium datasets, also extract score paths so users can switch between them
            dataset = VisiumDataset(
                dataset_id=dataset_id,
                alias=f"{dataset_type} Dataset",
                adata_path=adata_path,
                user=user,
                tangram_adata_path=output_files.get("adata_tg_scores_path") if use_tangram else None,
                geojson_path=geojson_path,
                genie_network_path=genie_network_path,
                sponge_network_path=sponge_network_path,
                created_at=datetime.now(),
                dataset_type=dataset_type,
                use_tangram=use_tangram,
                use_multiome=use_multiome,
                adata_st_scores_path=output_files.get("adata_st_scores_path"),
                adata_tg_scores_path=output_files.get("adata_tg_scores_path"),
                tf_graph_path=output_files.get("tf_graph_path")
            )

        # Register it
        self.register_uploaded_dataset(dataset=dataset)
        print(f"✓ Registered dataset from config: {dataset_id} (type: {type(dataset).__name__})")
        return dataset

    def _resolve_config_path(self, path_str: Optional[str], base_dir: Path) -> Optional[str]:
        """Resolve a path from config, whether absolute or relative to base_dir."""
        if not path_str:
            return None
        if path_str.startswith("/api/"):
            return path_str
        p = Path(path_str)
        if p.is_absolute() and p.exists():
            return str(p)
        candidate = base_dir / p
        if candidate.exists():
            return str(candidate)
        candidates = list(base_dir.rglob(p.name))
        if candidates:
            return str(candidates[0])
        return str(candidate if not p.is_absolute() else p)

    def register_builtin_from_config(
        self,
        config_file: Path,
        dataset_id: str,
        alias: str,
        description: Optional[str] = None,
    ) -> Dataset:
        """
        Register a built-in dataset directly from its configuration file.
        Infers all paths and scores relative to the dataset directory.
        """
        config_file = Path(config_file)
        if not config_file.exists():
            raise ValueError(f"Builtin config file not found: {config_file}")

        with open(config_file, 'r') as f:
            config = json.load(f)

        base_dir = config_file.parent
        output_files = config.get("output_files", {})
        dataset_type = config.get("dataset", "Visium")

        if dataset_type == "Visium" and config.get("multiome", {}).get("use") == True:
            dataset_type = "Multiome"

        # Determine adata path
        raw_adata_path = self._determine_adata_path(dataset_type, output_files, config)
        adata_path = self._resolve_config_path(raw_adata_path, base_dir)

        if not adata_path or not Path(adata_path).exists():
            # Fallback to finding any .h5ad
            h5ad_files = list(base_dir.rglob("*scores*.h5ad")) or list(base_dir.glob("*.h5ad"))
            if h5ad_files:
                adata_path = str(h5ad_files[0])
            else:
                raise ValueError(f"AnnData file not found for builtin dataset {dataset_id} in {base_dir}")

        geojson_path = f"/api/geojson/{dataset_id}"
        genie_network_path = self._resolve_config_path(output_files.get("genie_network_path"), base_dir)
        sponge_network_path = self._resolve_config_path(output_files.get("sponge_network_path"), base_dir)
        tangram_adata_path = self._resolve_config_path(output_files.get("adata_tg_scores_path"), base_dir)
        adata_st_scores_path = self._resolve_config_path(output_files.get("adata_st_scores_path"), base_dir)
        adata_tg_scores_path = self._resolve_config_path(output_files.get("adata_tg_scores_path"), base_dir)
        tf_graph_path = self._resolve_config_path(output_files.get("tf_graph_path"), base_dir)

        use_tangram = config.get("tangram", {}).get("use", False) if config.get("tangram") else False
        use_multiome = config.get("multiome", {}).get("use", False) if config.get("multiome") else False

        dataset_type_lower = dataset_type.lower()
        if "xenium" in dataset_type_lower:
            xenium_grid_adata_path = self._resolve_config_path(output_files.get("xenium_grid_adata_path") or raw_adata_path, base_dir)
            dataset = XeniumDataset(
                dataset_id=dataset_id,
                alias=alias,
                adata_path=adata_path,
                user="builtin",
                tangram_adata_path=tangram_adata_path if use_tangram else None,
                geojson_path=geojson_path,
                genie_network_path=genie_network_path,
                sponge_network_path=sponge_network_path,
                tf_graph_path=tf_graph_path,
                created_at=datetime.now(),
                dataset_type=dataset_type,
                use_tangram=use_tangram,
                use_multiome=use_multiome,
                xenium_grid_adata_path=xenium_grid_adata_path,
                description=description,
            )
        elif "multiome" in dataset_type_lower:
            dataset = MultiomeDataset(
                dataset_id=dataset_id,
                alias=alias,
                adata_path=adata_path,
                user="builtin",
                tangram_adata_path=adata_tg_scores_path,
                geojson_path=geojson_path,
                genie_network_path=genie_network_path,
                sponge_network_path=sponge_network_path,
                created_at=datetime.now(),
                adata_st_scores_path=adata_st_scores_path,
                adata_tg_scores_path=adata_tg_scores_path,
                adata_map_path=self._resolve_config_path(output_files.get("adata_map_path"), base_dir),
                adata_map_X_csv_path=self._resolve_config_path(output_files.get("adata_map_X_csv_path"), base_dir),
                adata_map_var_csv_path=self._resolve_config_path(output_files.get("adata_map_var_csv_path"), base_dir),
                calc_scores_log_path=self._resolve_config_path(output_files.get("calc_scores_log_path"), base_dir),
                global_motif_analysis_path=self._resolve_config_path(output_files.get("global_motif_analysis_path"), base_dir),
                motif_to_tf_csv_path=self._resolve_config_path(output_files.get("motif_to_tf_csv_path"), base_dir),
                spot_obj_chromvar_path=self._resolve_config_path(output_files.get("spot_obj_chromvar_path"), base_dir),
                spot_obj_footprints_path=self._resolve_config_path(output_files.get("spot_obj_footprints_path"), base_dir),
                dissociated_obj_footprints_path=self._resolve_config_path(output_files.get("dissociated_obj_footprints_path"), base_dir),
                chromvar_scores_csv_path=self._resolve_config_path(output_files.get("chromvar_scores_csv_path"), base_dir),
                diff_motif_activity_csv_paths=output_files.get("diff_motif_activity_csv_paths"),
                footprint_pdf_paths=output_files.get("footprint_pdf_paths"),
                tf_graph_path=tf_graph_path,
                grn_evaluation_name=output_files.get("grn_evaluation_name"),
                description=description,
            )
        else:
            dataset = VisiumDataset(
                dataset_id=dataset_id,
                alias=alias,
                adata_path=adata_path,
                user="builtin",
                tangram_adata_path=tangram_adata_path if use_tangram else None,
                geojson_path=geojson_path,
                genie_network_path=genie_network_path,
                sponge_network_path=sponge_network_path,
                created_at=datetime.now(),
                dataset_type=dataset_type,
                use_tangram=use_tangram,
                use_multiome=use_multiome,
                adata_st_scores_path=adata_st_scores_path,
                adata_tg_scores_path=adata_tg_scores_path,
                tf_graph_path=tf_graph_path,
                description=description,
            )

        self.register_builtin_dataset(dataset=dataset)
        print(f"✓ Registered builtin dataset from config: {dataset_id} ({alias})")
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
                                tf_graph_path = config["output_files"].get("tf_graph_path")
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
                    tf_graph_path=tf_graph_path,
                )

                # Register the Dataset object
                self.register_uploaded_dataset(dataset=dataset)
                results[dataset_id] = "registered"

            except Exception as e:
                results[str(h5ad_file)] = f"error: {e}"
                import traceback
                traceback.print_exc()

        return results
