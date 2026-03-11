# like Diana's map_cells.py file, but slight structure changes

#!/usr/bin/env python3
# coding: utf-8

import argparse
import os
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.stats import spearmanr
import tangram as tg
import gene_selection
import torch
import math
# Allow writing nullable strings in anndata
import anndata as ad
ad.settings.allow_write_nullable_strings = True

def select_genes(ad_sc, ad_sp, selection_mode: str, cell_label: str):
    """
    Select genes according to the requested strategy. Falls back to overlapping genes
    when no selection is provided.
    """
    if selection_mode is not None:
        selection_mode = selection_mode.lower()

    # rank_genes_groups (used by ctg/spapros) requires the groupby column to be
    # categorical; cast it here so the caller doesn't have to worry about dtype.
    if cell_label and cell_label in ad_sc.obs.columns:
        ad_sc.obs[cell_label] = ad_sc.obs[cell_label].astype("category")

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


def run_tangram(ad_sc: object, ad_sp: object, gene_selection_mode: str = None, cell_label: str = 'cell_type', ensembl_col: str = "", feature_col: str = "", device_choice: str = 'cpu', multiome: bool = False):
    # Tangram preprocessing mutates ad_sc/ad_sp in-place.
    # We keep copies only to restore spatial metadata + map var annotations later.
    adata_sp_copy = ad_sp.copy()
    adata_sc_copy = ad_sc.copy()
    ad_sc_use = ad_sc.raw.to_adata().copy() if ad_sc.raw is not None else ad_sc.copy()
    # if multiome -> mode = cells
    if multiome:
        mode = "cells"
        print("[mapping] Multiome mode: mapping individual cells to space (mode='cells').")
    else:
        mode = "clusters"
        print("[mapping] Non-multiome mode: mapping clusters to space (mode='clusters').")


    # Gene selection
    genes = select_genes(ad_sc_use, ad_sp, gene_selection_mode, cell_label)

    # Preprocessing
    # split into train and tet to evaluate the mapping performance on the test set
    tg.pp_adatas(ad_sc_use, ad_sp, genes=genes)

    # Mapping logic based on user choice
    if device_choice == "gpu":
        if not torch.cuda.is_available():
            #raise RuntimeError("GPU selected but no CUDA device is available.")
            print("[mapping] Using CPU, no CUDA device available")
            ad_map = tg.map_cells_to_space(ad_sc_use, ad_sp, mode=mode, device="cpu", cluster_label=cell_label)
            if mode == "clusters":
                ad_ge = tg.project_genes(ad_map, ad_sc_use, cluster_label=cell_label)
            else:
                ad_ge = tg.project_genes(ad_map, ad_sc_use)
        else:
            print("[mapping] Using GPU (cuda:0) with cluster_label.")
            # M = ad_map.X map dissociated to spatial
            ad_map = tg.map_cells_to_space(ad_sc_use, ad_sp, mode=mode, device="cuda:0", cluster_label=cell_label)
            if mode == "clusters":
                ad_ge = tg.project_genes(ad_map, ad_sc_use, cluster_label=cell_label)
            else:
                ad_ge = tg.project_genes(ad_map, ad_sc_use)
    elif device_choice == "cpu":
        print("[mapping] Using CPU with mode='clusters' and no cluster_label.")
        ad_map = tg.map_cells_to_space(ad_sc_use, ad_sp, mode=mode, device="cpu", cluster_label=cell_label)
        if mode == "clusters":
            ad_ge = tg.project_genes(ad_map, ad_sc_use, cluster_label=cell_label)
        else:
            ad_ge = tg.project_genes(ad_map, ad_sc_use)
    else:
        raise ValueError("Invalid device choice. Use 'cpu' or 'gpu'.")

    # Project gene expression and cell annotations
    print("[mapping] Projecting gene expression and cell annotations to space...")
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

def spotwise_gene_correlations(ad_ge_pred, ad_sp_true, genes):
    common = [g for g in genes if g in ad_ge_pred.var_names and g in ad_sp_true.var_names]
    if len(common) == 0:
        raise ValueError("No common genes found for evaluation.")

    Xp = ad_ge_pred[:, common].X
    Xt = ad_sp_true[:, common].X

    Xp = Xp.A if hasattr(Xp, "A") else np.asarray(Xp)
    Xt = Xt.A if hasattr(Xt, "A") else np.asarray(Xt)

    corrs = []
    for j in range(len(common)):
        # correlation across spots for gene j
        c = spearmanr(Xp[:, j], Xt[:, j]).correlation
        corrs.append(c)

    corrs = np.array(corrs, dtype=float)
    return common, corrs


def mapping_entropy_and_mass(ad_map, eps=1e-12):
    """
     Confidence metrics for the mapping:
     Calculate the entropy of the mapping for each row (cell or cluster) and the total mass assigned to each spot.
     - Very high entropy for most rows = mapping is uncertain / diffuse (often batch mismatch).
     - Extremely uneven spot_mass = collapse (model dumps everything into few spots).
     """

    M = ad_map.X
    M = M.A if hasattr(M, "A") else np.asarray(M)

    # already row-normalized M i think
    # row_sums = M.sum(axis=1, keepdims=True) + eps
    # M = M / row_sums

    entropy_per_row = -np.sum(M * np.log(M + eps), axis=1)  # rows = cells or clusters
    spot_mass = M.sum(axis=0)  # total mass per spot

    return entropy_per_row, spot_mass

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sc_path", required=True, help="Pfad zur Single-Cell h5ad")
    parser.add_argument("--sp_path", required=True, help="Pfad zur Spatial h5ad")
    parser.add_argument("--multiome", action='store_true', help='Whether the single-cell data is from a multiome experiment (default: False)')
    parser.add_argument("--outdir", required=True, help="Pfad für Ausgabe-Verzeichnis (Tangram)")
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
    ad_ge, adata_map = run_tangram(
        ad_sc,
        ad_sp,
        multiome=args.multiome,
        gene_selection_mode=args.gene_selection_mode,
        cell_label=args.cell_label,
        ensembl_col=args.ensembl_col,
        feature_col=args.feature_col,
        device_choice=args.device,
    )

    # Ergebnis speichern
    print("[tangram_cli] Schreibe Ergebnis nach:", args.outdir)
    print(ad_ge)
    print(ad_ge.obsm["tangram_ct_pred"])
    ad_ge.obsm["tangram_ct_pred"] = ad_ge.obsm["tangram_ct_pred"].rename(columns=str)
    ad_ge.write_h5ad(os.path.join(args.outdir, "tangram_output.h5ad"))
    adata_map.write(os.path.join(args.outdir, "adata_map.h5ad"))
    # wtite also X as csv and var/obs as csv

    df_adata_map = pd.DataFrame(
        adata_map.X,  # convert sparse to dense
        index=adata_map.obs_names,
        columns=adata_map.var_names)
    df_adata_map.to_csv(os.path.join(args.outdir, "adata_map.X.csv"))
    #mmwrite(os.path.join(args.outdir, "adata_map.X.mtx"), adata_map.X)
    adata_map.var.to_csv(os.path.join(args.outdir, "adata_map.var.csv"))
    adata_map.obs.to_csv(os.path.join(args.outdir, "adata_map.obs.csv"))
    print("[tangram_cli] Fertig.")


