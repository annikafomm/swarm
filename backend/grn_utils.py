"""
GRN Evaluation Utilities
Handles file I/O, graph filtering, CSV parsing, and directory scanning for GRN evaluation results.
"""

from pathlib import Path
from typing import Any, Dict, List, Optional
import json
import csv
import os


class GrnEvaluationUtils:
    """Utilities for GRN Evaluation file I/O and graph filtering"""

    @staticmethod
    def list_available_runs(job_dir: Path) -> List[str]:
        """
        List all available on-demand GRN evaluation runs by scanning the on_demand directory.
        Returns sorted list of run names (directory names).
        """
        on_demand_dir = job_dir / "multiome" / "GRN_evaluation" / "on_demand"
        if not on_demand_dir.exists():
            return []

        try:
            runs = [d.name for d in on_demand_dir.iterdir() if d.is_dir()]
            return sorted(runs)
        except Exception as e:
            print(f"Error listing GRN runs: {e}")
            return []

    @staticmethod
    def load_graph_json(json_path: Path) -> Optional[Dict[str, Any]]:
        """Load graph.json file"""
        if not json_path.exists():
            return None

        try:
            with open(json_path, 'r') as f:
                return json.load(f)
        except Exception as e:
            print(f"Error loading graph JSON {json_path}: {e}")
            return None

    @staticmethod
    def filter_graph_by_mode(
        graph_data: Dict[str, Any],
        filter_mode: str = "full"
    ) -> Dict[str, Any]:
        """
        Filter graph based on mode:
        - 'prior': Keep only edges from priorTF nodes
        - 'extended': Remove red/yellow edges (+ orphaned nodes)
        - 'full': No filtering (default)

        Returns filtered copy of graph_data with nodes and links arrays
        """
        if not graph_data or "nodes" not in graph_data or "links" not in graph_data:
            return graph_data

        # Create shallow copy to avoid modifying original
        filtered = {
            "nodes": [n.copy() for n in graph_data.get("nodes", [])],
            "links": [l.copy() for l in graph_data.get("links", [])]
        }

        if filter_mode == "full":
            return filtered

        # Get all node IDs
        all_node_ids = {node["id"] for node in filtered["nodes"]}

        if filter_mode == "prior":
            # Keep only edges from priorTF nodes
            prior_tf_ids = {
                node["id"] for node in filtered["nodes"]
                if node.get("priorTF") is True
            }
            filtered["links"] = [
                link for link in filtered["links"]
                if link.get("source") in prior_tf_ids
            ]

        elif filter_mode == "extended":
            # Remove red and yellow edges
            filtered["links"] = [
                link for link in filtered["links"]
                if link.get("edge_color") not in ["red", "yellow"]
            ]

        # Find nodes referenced in filtered links
        referenced_node_ids = set()
        for link in filtered["links"]:
            src = link.get("source")
            tgt = link.get("target")
            if src is not None:
                referenced_node_ids.add(src)
            if tgt is not None:
                referenced_node_ids.add(tgt)

        # Keep target nodes (type="to") and only referenced nodes
        target_node_ids = {
            node["id"] for node in filtered["nodes"]
            if node.get("type") == "to"
        }

        filtered["nodes"] = [
            node for node in filtered["nodes"]
            if node["id"] in referenced_node_ids or node["id"] in target_node_ids
        ]

        return filtered

    @staticmethod
    def list_plot_images(plots_dir: Path) -> Dict[str, List[str]]:
        """
        List PNG images organized by subdirectory.
        Returns: {"peak_plots": [...], "motif_plots": [...]}
        """
        result = {}
        if not plots_dir.exists():
            return result

        try:
            for subdir in ["peak_plots", "motif_plots"]:
                subdir_path = plots_dir / subdir
                if subdir_path.exists():
                    images = sorted([
                        img_file.name for img_file in subdir_path.glob("*.png")
                    ])
                    if images:
                        result[subdir] = images
        except Exception as e:
            print(f"Error listing plots: {e}")

        return result

    @staticmethod
    def read_csv_as_json(csv_path: Path) -> Optional[Dict[str, Any]]:
        """
        Read CSV file and return as JSON-serializable format.
        Attempts to convert numeric columns to floats.
        Returns: {"rows": [...], "columns": [...]}
        """
        if not csv_path.exists():
            return None

        try:
            rows = []
            columns = []

            with open(csv_path, 'r', encoding='utf-8') as f:
                reader = csv.DictReader(f)
                columns = reader.fieldnames or []

                for row in reader:
                    converted_row = {}
                    for key, value in row.items():
                        # Try to convert to float, otherwise keep as string
                        try:
                            converted_row[key] = float(value)
                        except (ValueError, TypeError):
                            converted_row[key] = value
                    rows.append(converted_row)

            return {
                "rows": rows,
                "columns": list(columns)
            }
        except Exception as e:
            print(f"Error reading CSV {csv_path}: {e}")
            return None

    @staticmethod
    def get_grn_evaluation_path(
        job_dir: Path,
        evaluation_name: str,
        source: str = "both"
    ) -> Optional[Path]:
        """
        Get the path to a GRN evaluation directory.
        source can be: 'initial_upload', 'on_demand', or 'both' (tries both)
        """
        if source in ["on_demand", "both"]:
            path = job_dir / "multiome" / "GRN_evaluation" / "on_demand" / evaluation_name
            if path.exists():
                return path

        if source in ["initial_upload", "both"]:
            path = job_dir / "multiome" / "GRN_evaluation" / "initial_upload" / evaluation_name
            if path.exists():
                return path

        return None

    @staticmethod
    def get_grn_stats_from_meta(meta: Dict[str, Any], eval_type: str, eval_name: str) -> Optional[Dict[str, Any]]:
        """
        Extract stats from meta dict for initial_upload scenario.
        eval_type: 'peak' or 'motif'
        eval_name: GRN evaluation name
        Returns stats as table format or None
        """
        type_key = f"{eval_type}_stats"
        if type_key not in meta:
            return None

        stats_by_name = meta[type_key]
        if not isinstance(stats_by_name, dict):
            return None

        stats = stats_by_name.get(eval_name)
        if not stats:
            return None

        # Convert to standard format
        if isinstance(stats, dict):
            # Assume it's already in {rows: [...], columns: [...]} format or similar
            return stats

        return None
