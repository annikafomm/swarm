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
    if not selection_mode is None:
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

    if not genes is None:
        genes = list(genes)
        print(f"[gene_selection] Selected n={len(genes)} genes")
    else:
        print(f"[gene_selection] Selected n=0 genes")
    return genes


def run_tangram(ad_sc: object, ad_sp: object, gene_selection_mode: str = None, cell_label: str = 'cell_type', ensembl_col: str = "", feature_col: str = "", device_choice: str = 'cpu'):
    adata_sp_copy = ad_sp.copy()
    adata_sc_copy = ad_sc.copy()

    # Gene selection
    genes = select_genes(ad_sc, ad_sp, gene_selection_mode, cell_label)

    # Preprocessing
    tg.pp_adatas(ad_sc, ad_sp, genes=genes)

    # Mapping logic based on user choice
    if device_choice == "gpu":
        if not torch.cuda.is_available():
            #raise RuntimeError("GPU selected but no CUDA device is available.")
            print("[mapping] Using CPU, no CUDA device available")
            ad_map = tg.map_cells_to_space(ad_sc, ad_sp, device="cpu")
            ad_ge = tg.project_genes(ad_map, ad_sc)
        else:
            print("[mapping] Using GPU (cuda:0) with cluster_label.")
            # M = ad_map.X map dissociated to spatial
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

    return(ad_ge, ad_map)

