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

    def __init__(self, registry_file: Path = Path("./uploads/dataset_registry.json")):
        self.registry_file = registry_file
        self.registry_file.parent.mkdir(parents=True, exist_ok=True)
        self.datasets = self._load_registry()

    def _load_registry(self) -> Dict:
        """Load existing registry or create new one"""
        if self.registry_file.exists():
            try:
                with open(self.registry_file, 'r') as f:
                    return json.load(f)
            except Exception as e:
                print(f"Error loading registry: {e}")
                return {"builtin": {}, "uploaded": {}}
        return {"builtin": {}, "uploaded": {}}

    def _save_registry(self):
        """Persist registry to disk"""
        with open(self.registry_file, 'w') as f:
            json.dump(self.datasets, indent=2, fp=f)

    def register_uploaded_dataset(
        self,
        dataset_id: str,
        alias: str,
        adata_path: str,
        tangram_adata_path: Optional[str] = None,
        genie_network_path: Optional[str] = None,
        sponge_network_path: Optional[str] = None,
        user: str = "anonymous",
        **metadata
    ) -> None:
        """Register a newly uploaded dataset"""

        if "uploaded" not in self.datasets:
            self.datasets["uploaded"] = {}

        self.datasets["uploaded"][dataset_id] = {
            "alias": alias,
            "adata_path": adata_path,
            "tangram_adata_path": tangram_adata_path,
            "genie_network_path": genie_network_path,
            "sponge_network_path": sponge_network_path,
            "user": user,
            "created_at": datetime.now().isoformat(),
            **metadata
        }
        self._save_registry()

    def register_builtin_dataset(
        self,
        dataset_id: str,
        alias: str,
        adata_path: str,
        **metadata
    ) -> None:
        """Register a built-in dataset (from backend/data)"""

        if "builtin" not in self.datasets:
            self.datasets["builtin"] = {}

        self.datasets["builtin"][dataset_id] = {
            "alias": alias,
            "adata_path": adata_path,
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
