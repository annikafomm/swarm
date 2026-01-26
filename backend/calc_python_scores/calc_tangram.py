# like Diana's map_cells.py file, but slight structure changes

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
import math


def select_genes(ad_sc, ad_sp, selection_mode: str, cell_label: str):
    """
    Select genes according to the requested strategy. Falls back to overlapping genes
    when no selection is provided.
    """
    if selection_mode is not None:
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
        genes = None
        print("[gene_selection] Using all overlapping genes (no selection mode)")

    if genes is None:
        return None

    genes = list(genes)
    print(f"[gene_selection] Selected n={len(genes)} genes")
    return genes


def run_tangram(ad_sc: object, ad_sp: object, gene_selection_mode: str = None, cell_label: str = 'cell_type', ensembl_col: str = "", feature_col: str = "", device_choice: str = 'cpu'):
    # Tangram preprocessing mutates ad_sc/ad_sp in-place.
    # We keep copies only to restore spatial metadata + map var annotations later.
    adata_sp_copy = ad_sp.copy()
    adata_sc_copy = ad_sc.copy()

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
        print("[mapping] Using CPU with mode='clusters' and with cluster_label.")
        ad_map = tg.map_cells_to_space(ad_sc, ad_sp, mode="clusters", device="cpu", cluster_label=cell_label, num_epochs=200)
        ad_ge = tg.project_genes(ad_map, ad_sc, cluster_label=cell_label)
    else:
        raise ValueError("Invalid device choice. Use 'cpu' or 'gpu'.")

    # Project gene expression and cell annotations
    tg.project_cell_annotations(ad_map, ad_ge, annotation=cell_label)

    ad_ge.obsm['spatial'] = adata_sp_copy.obsm['spatial']
    ad_ge.uns['spatial'] = adata_sp_copy.uns['spatial']

    if (ensembl_col != "" and feature_col != "") and (ensembl_col in adata_sc_copy.var.columns and feature_col in adata_sc_copy.var.columns):

        # Create mapping dictionaries with uppercase gene names for case-insensitive matching
        ensembl_map = dict(zip([x.upper() for x in adata_sc_copy.var_names], adata_sc_copy.var[ensembl_col]))
        feature_map = dict(zip([x.upper() for x in adata_sc_copy.var_names], adata_sc_copy.var[feature_col]))

        # Assign columns using case-insensitive mapping
        ad_ge.var[ensembl_col] = [ensembl_map.get(x.upper(), pd.NA) for x in ad_ge.var_names]
        ad_ge.var[feature_col] = [feature_map.get(x.upper(), pd.NA) for x in ad_ge.var_names]


    print(f"[done] using label '{cell_label}'.")

    return(ad_ge)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sc_path", required=True, help="Pfad zur Single-Cell h5ad")
    parser.add_argument("--sp_path", required=True, help="Pfad zur Spatial h5ad")
    parser.add_argument("--out_path", required=True, help="Pfad für Ausgabe-h5ad (Tangram)")
    parser.add_argument("--gene_selection_mode", default=None, help="ctg/hvg/spapros/svg/None")
    parser.add_argument("--cell_label", default="cell_type")
    parser.add_argument("--ensembl_col", default="")
    parser.add_argument("--feature_col", default="")
    parser.add_argument("--device", default="cpu", choices=["cpu", "gpu"])

    args = parser.parse_args()

    # Daten laden
    print("[tangram_cli] Lade Single-Cell AnnData:", args.sc_path)
    ad_sc = sc.read_h5ad(args.sc_path)

    print("[tangram_cli] Lade Spatial AnnData:", args.sp_path)
    ad_sp = sc.read_h5ad(args.sp_path)

    # Tangram ausführen
    ad_ge = run_tangram(
        ad_sc,
        ad_sp,
        gene_selection_mode=args.gene_selection_mode,
        cell_label=args.cell_label,
        ensembl_col=args.ensembl_col,
        feature_col=args.feature_col,
        device_choice=args.device,
    )

    # Ergebnis speichern
    print("[tangram_cli] Schreibe Ergebnis nach:", args.out_path)
    ad_ge.write_h5ad(args.out_path)
    print("[tangram_cli] Fertig.")


