#!/usr/bin/env python3
# coding: utf-8

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
import tangram as tg
import gene_selection
import torch


def load_data(sc_path: str, sp_path: str):
    """Load single-cell and spatial AnnData objects."""
    ad_sp = sc.read_h5ad(sp_path)
    ad_sc = sc.read_h5ad(sc_path)
    return ad_sc, ad_sp


def select_genes(ad_sc, ad_sp, selection_mode: str, cell_label: str):
    """
    Select genes according to the requested strategy. Falls back to overlapping genes
    when no selection is provided.
    """
    selection_mode = selection_mode.lower()
    if selection_mode == "ctg":
        genes = gene_selection.ctg(ad_sc, cell_label)
        print("[gene_selection] ctg")
    elif selection_mode == "hvg":
        genes = gene_selection.hvg(ad_sc)
        print("[gene_selection] hvg")
    elif selection_mode == "spapros":
        sc.pp.highly_variable_genes(ad_sc, flavor="seurat", n_top_genes=2000)
        genes = gene_selection.spapros(ad_sc, cell_label)
        print("[gene_selection] spapros")
    elif selection_mode == "svg":
        ad_sp.raw = ad_sp.copy()
        genes = gene_selection.svg(ad_sp)
        print("[gene_selection] svg")
    else:
        print("No gene selection performed. All overlapping genes will be used.")
        genes = None
        print("none")

    genes = list(genes)
    print(f"[gene_selection] Selected n={len(genes)} genes")
    return genes


def ensure_outdir(path: str):
    os.makedirs(path, exist_ok=True)


def run_pipeline(sc_path: str, sp_path: str, output_dir: str, gene_selection_mode: str, cell_label: str, device_choice: str):
    # Load
    ad_sc, ad_sp = load_data(sc_path, sp_path)

    # Gene selection
    genes = select_genes(ad_sc, ad_sp, gene_selection_mode, cell_label)

    # Preprocessing
    tg.pp_adatas(ad_sc, ad_sp, genes=genes)

    # Mapping logic based on user choice
    if device_choice == "gpu":
        if not torch.cuda.is_available():
            raise RuntimeError("GPU selected but no CUDA device is available.")
        print("[mapping] Using GPU (cuda:0) with cluster_label.")
        ad_map = tg.map_cells_to_space(ad_sc, ad_sp, device="cuda:0")
        ad_ge = tg.project_genes(ad_map, ad_sc)
    elif device_choice == "cpu":
        print("[mapping] Using CPU with mode='clusters' and no cluster_label.")
        ad_map = tg.map_cells_to_space(ad_sc, ad_sp, mode="clusters", device="cpu", cluster_label=cell_label)
        ad_ge = tg.project_genes(ad_map, ad_sc, cluster_label=cell_label)
    else:
        raise ValueError("Invalid device choice. Use 'cpu' or 'gpu'.")

    # Project gene expression and cell annotations
    tg.project_cell_annotations(ad_map, ad_ge, annotation=cell_label)

    # Write outputs
    ensure_outdir(output_dir)
    ad_map.write(os.path.join(output_dir, f"ad_map_{cell_label}.h5ad"))
    ad_ge.write(os.path.join(output_dir, f"ad_ge_{cell_label}.h5ad"))
    print(f"[done] Wrote ad_map and ad_ge to '{output_dir}' using label '{cell_label}'.")


def build_argparser():
    p = argparse.ArgumentParser(description="Merged Tangram cell mapping pipeline (manual CPU/GPU selection).")
    p.add_argument('--sc_path', type=str, required=True, help="Path to the single-cell .h5ad file.")
    p.add_argument('--sp_path', type=str, required=True, help="Path to the spatial .h5ad file.")
    p.add_argument('--output_dir', type=str, required=True, help="Directory for outputs.")
    p.add_argument('--gene_selection', type=str,
                   choices=['ctg', 'hvg', 'spapros', 'svg'],
                   default=None,
                   help="Gene selection strategy. Default: use all overlapping genes.")
    p.add_argument('--cell_label', type=str, default='cell_type',
                   help="Column in ad_sc.obs with cluster/cell annotations (e.g. 'cell_type' or 'cell_subclass').")
    p.add_argument('--device', type=str, choices=['cpu', 'gpu'], required=True,
                   help="Select whether to run on CPU or GPU.")
    return p


def main():
    args = build_argparser().parse_args()
    run_pipeline(args.sc_path, args.sp_path, args.output_dir, args.gene_selection, args.cell_label, args.device)


if __name__ == "__main__":
    main()
