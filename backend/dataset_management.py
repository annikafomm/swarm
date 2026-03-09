import json
import os
from pathlib import Path
from datetime import datetime
from typing import Optional, Dict, List

class DatasetRegistry:
    """
    Manages available datasets across the application.
    Stores metadata about uploaded and built-in datasets.
    """

    def __init__(self, registry_file: Path = None):
        if registry_file is None:
            # Use absolute path based on script location
            base_path = Path(__file__).parent
            registry_file = base_path / "uploads" / "dataset_registry.json"
        self.registry_file = registry_file
        self.registry_file.parent.mkdir(parents=True, exist_ok=True)
        print(f"Dataset registry path: {self.registry_file}")
        self.datasets = self._load_registry()

    def _load_registry(self, clear_uploads: bool = True) -> Dict:
        """Load existing registry or create new one.

        Args:
            clear_uploads: If True (on startup), clear uploaded datasets.
                          If False (runtime), preserve existing uploads.
        """
        if self.registry_file.exists():
            try:
                with open(self.registry_file, 'r') as f:
                    data = json.load(f)
                    if clear_uploads:
                        # On startup, keep only builtin datasets, clear uploaded ones
                        return {
                            "builtin": data.get("builtin", {}),
                            "uploaded": {}
                        }
                    else:
                        # At runtime, preserve everything
                        return data

            except Exception as e:
                print(f"Error loading registry: {e}")
                return {"builtin": {}, "uploaded": {}}
        return {"builtin": {}, "uploaded": {}}

    def _save_registry(self):
        """Persist registry to disk"""
        try:
            with open(self.registry_file, 'w') as f:
                json.dump(self.datasets, indent=2, fp=f)
            print(f"✓ Registry saved to {self.registry_file}")
        except Exception as e:
            print(f"✗ Error saving registry: {e}")
            raise

    def register_uploaded_dataset(
        self,
        dataset_id: str,
        alias: str,
        adata_path: str,
        tangram_adata_path: Optional[str] = None,
        geojson_path: Optional[str] = None,
        genie_network_path: Optional[str] = None,
        sponge_network_path: Optional[str] = None,
        user: str = "anonymous",
        footprint_list: Optional[List[str]] = None,
        **metadata
    ) -> None:
        """Register a newly uploaded dataset"""

        if "uploaded" not in self.datasets:
            self.datasets["uploaded"] = {}

        self.datasets["uploaded"][dataset_id] = {
            "alias": alias,
            "adata_path": adata_path,
            "tangram_adata_path": tangram_adata_path,
            "geojson_path": geojson_path,
            "genie_network_path": genie_network_path,
            "sponge_network_path": sponge_network_path,
            "user": user,
            "created_at": datetime.now().isoformat(),
            "footprint_list": footprint_list,
        }
        print(f"Saving registry to {self.registry_file}")
        self._save_registry()
        print(f"✓ Registered dataset {dataset_id}")

    def register_builtin_dataset(
        self,
        dataset_id: str,
        alias: str,
        adata_path: str,
        geojson_path: Optional[str] = None,
        genie_network_path: Optional[str] = None,
        sponge_network_path: Optional[str] = None,
        **metadata
    ) -> None:
        """Register a built-in dataset (from backend/data)"""

        if "builtin" not in self.datasets:
            self.datasets["builtin"] = {}

        self.datasets["builtin"][dataset_id] = {
            "alias": alias,
            "adata_path": adata_path,
            "geojson_path": geojson_path,
            "genie_network_path": genie_network_path,
            "sponge_network_path": sponge_network_path,
            "created_at": datetime.now().isoformat(),
            **metadata
        }
        self._save_registry()

    def get_all_datasets(self) -> Dict[str, Dict]:
        """Return all available datasets grouped by type"""
        return self.datasets

    def get_dataset_by_id(self, dataset_id: str) -> Optional[Dict]:
        """Get a single dataset by ID"""
        for category in ["builtin", "uploaded"]:
            if category in self.datasets and dataset_id in self.datasets[category]:
                return self.datasets[category][dataset_id]
        return None

    def get_user_datasets(self, user: str) -> Dict[str, Dict]:
        """Get all datasets for a specific user"""
        user_datasets = {}
        if "uploaded" in self.datasets:
            for dataset_id, info in self.datasets["uploaded"].items():
                if info.get("user") == user:
                    user_datasets[dataset_id] = info
        return user_datasets

    def remove_dataset(self, dataset_id: str) -> bool:
        """Remove a dataset from registry"""
        for category in ["builtin", "uploaded"]:
            if category in self.datasets and dataset_id in self.datasets[category]:
                del self.datasets[category][dataset_id]
                self._save_registry()
                return True
        return False

    def update_dataset_paths(
        self,
        dataset_id: str,
        **new_paths
    ) -> bool:
        """Update paths for an existing dataset"""
        dataset = self.get_dataset_by_id(dataset_id)
        if dataset:
            dataset.update(new_paths)
            self._save_registry()
            return True
        return False

    def rescan_uploads_folder(self, uploads_dir: Path) -> Dict[str, str]:
        """
        TESTING ONLY: Scan the uploads folder and re-register datasets found there.
        Looks for *_st_scores.h5ad files and attempts to reconstruct metadata.
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

                if config_file.exists():
                    try:
                        with open(config_file, 'r') as f:
                            config = json.load(f)
                            metadata = {
                                "dataset_type": config.get("dataset", "Unknown"),
                                "use_tangram": config.get("tangram", {}).get("use", False),
                            }
                            if "output_files" in config:
                                genie_network_path = config["output_files"].get("genie_network_path")
                                sponge_network_path = config["output_files"].get("sponge_network_path")
                    except Exception as e:
                        print(f"  ⚠ Could not parse config for {dataset_id}: {e}")

                # Look for geojson in multiple locations
                # First in h5ad's direct directory, then in parent directories up to job_dir
                geojson_candidates = []
                search_dir = h5ad_file.parent
                max_depth = 5  # Limit how far up we search
                for _ in range(max_depth):
                    candidates = list(search_dir.glob("hexagons.geojson"))
                    if candidates:
                        geojson_candidates = candidates
                        geojson_path = f"/api/geojson/{dataset_id}"
                        break
                    # Stop if we've reached or passed the job_dir
                    if search_dir <= job_dir:
                        break
                    # Move up one directory level
                    search_dir = search_dir.parent

                # Extract user from directory name (job_TIMESTAMP_USER format)
                user_part = dataset_id.split("_", 2)[-1] if "_" in dataset_id else "unknown"

                footprint_list = [
                    str(h5ad_file.parent.relative_to(uploads_dir) / f)
                    for f in os.listdir(h5ad_file.parent)
                    if f.startswith("footprint") and f.endswith(".pdf")
                ]
                # Register the dataset
                self.register_uploaded_dataset(
                    dataset_id=dataset_id,
                    alias=f"[RESCANNED] {dataset_id}",
                    adata_path=str(h5ad_file),
                    geojson_path=geojson_path,
                    genie_network_path=genie_network_path,
                    sponge_network_path=sponge_network_path,
                    user=user_part,
                    footprint_list=footprint_list,
                    **metadata
                )
                results[dataset_id] = "registered"

            except Exception as e:
                results[str(h5ad_file)] = f"error: {e}"

        return results
