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
import spapros as sp
from scipy import sparse
from SpatialDE.anndata import spatialde_test
import re
# Allow writing nullable strings in anndata
import anndata as ad
ad.settings.allow_write_nullable_strings = True

def parse_selection_modes(selection_mode):
    """
    Accepts:
      - None
      - "ctg"
      - "ctg,hvg,svg"
      - ["ctg", "hvg"]
      - ["ctg,hvg", "svg"]
    Returns a unique list while preserving order.
    """
    if selection_mode is None:
        return []

    if isinstance(selection_mode, str):
        raw = [selection_mode]
    else:
        raw = list(selection_mode)

    modes = []
    for item in raw:
        if item is None:
            continue
        for part in str(item).split(","):
            part = part.strip().lower()
            if part and part != "none":
                modes.append(part)

    return list(dict.fromkeys(modes))


def load_user_genes(path: str, gene_column: str = ""):
    """
    Load user-provided genes from:
      - txt: one gene per line
      - csv: first column or named column
      - tsv: first column or named column
    """
    if not path:
        return []

    if not os.path.exists(path):
        raise FileNotFoundError(f"Gene list file not found: {path}")

    lower = path.lower()

    if lower.endswith(".csv"):
        df = pd.read_csv(path)
        if gene_column:
            if gene_column not in df.columns:
                raise ValueError(
                    f"Column '{gene_column}' not found in {path}. "
                    f"Available columns: {list(df.columns)}"
                )
            genes = df[gene_column].tolist()
        else:
            genes = df.iloc[:, 0].tolist()

    elif lower.endswith(".tsv"):
        df = pd.read_csv(path, sep="\t")
        if gene_column:
            if gene_column not in df.columns:
                raise ValueError(
                    f"Column '{gene_column}' not found in {path}. "
                    f"Available columns: {list(df.columns)}"
                )
            genes = df[gene_column].tolist()
        else:
            genes = df.iloc[:, 0].tolist()

    else:
        with open(path, "r", encoding="utf-8") as f:
            genes = [
                line.strip()
                for line in f
                if line.strip() and not line.strip().startswith("#")
            ]

    genes = [str(g).strip() for g in genes if pd.notna(g) and str(g).strip()]
    genes = list(dict.fromkeys(genes))

    print(f"[gene_selection] Loaded n={len(genes)} user genes from {path}")
    return genes

def select_genes(ad_sc, ad_sp, selection_mode, cell_label: str, user_genes=None):
    """
    Build the union of genes from:
      - one or more selection methods
      - optional user-provided gene list

    If nothing is provided, returns None so Tangram uses all overlapping genes.
    """
    modes = parse_selection_modes(selection_mode)
    user_genes = user_genes or []

    if cell_label and cell_label in ad_sc.obs.columns:
        ad_sc.obs[cell_label] = ad_sc.obs[cell_label].astype("category")

    selected_sources = {}

    for mode in modes:
        if mode == "ctg":
            genes = gene_selection.ctg(ad_sc, cell_label)
            print("[gene_selection] ctg")
        elif mode == "hvg":
            genes = gene_selection.hvg(ad_sc)
            print("[gene_selection] hvg")
        elif mode == "spapros":
            sc.pp.highly_variable_genes(ad_sc, flavor="seurat", n_top_genes=2000)

            selector = sp.se.ProbesetSelector(
                ad_sc,
                celltype_key=cell_label,
            )
            selector.select_probeset()

            probeset = selector.probeset
            if "selection" in probeset.columns:
                genes = probeset.index[probeset["selection"].astype(bool)].tolist()
            else:
                genes = probeset.index.tolist()

            print("[gene_selection] spapros")
        elif mode == "svg":
            # NOT TESTED YET
            if ad_sp.raw is None and "counts" not in ad_sp.layers:
                print("[gene_selection] svg skipped: no raw counts available for SpatialDE")
                genes = []
            else:
                # Prefer raw counts from layers["counts"], otherwise fall back to ad_sp.raw
                if "counts" in ad_sp.layers:
                    ad_sp_tmp = ad_sp.copy()
                    ad_sp_tmp.X = ad_sp.layers["counts"].copy()
                else:
                    ad_sp_tmp = ad_sp.raw.to_adata().copy()
                    if "spatial" in ad_sp.obsm:
                        ad_sp_tmp.obsm["spatial"] = ad_sp.obsm["spatial"]

                if "spatial" not in ad_sp_tmp.obsm:
                    raise ValueError("ad_sp.obsm['spatial'] is required for SpatialDE.")

                # spatialde_test / NaiveDE expects a dense matrix here
                if sparse.issparse(ad_sp_tmp.X):
                    ad_sp_tmp.X = ad_sp_tmp.X.toarray()

                # Add spatial coordinates and library size covariate
                ad_sp_tmp.obs["x"] = np.asarray(ad_sp_tmp.obsm["spatial"][:, 0]).ravel()
                ad_sp_tmp.obs["y"] = np.asarray(ad_sp_tmp.obsm["spatial"][:, 1]).ravel()
                ad_sp_tmp.obs["total_counts"] = np.asarray(ad_sp_tmp.X.sum(axis=1)).ravel()

                # Run SpatialDE on raw counts
                svg_full = spatialde_test(
                    ad_sp_tmp,
                    coord_columns=["x", "y"],
                    regress_formula="~np.log(total_counts)",
                )

                # Extract selected genes from the result table
                if "qval" in svg_full.columns and "g" in svg_full.columns:
                    genes = svg_full.loc[svg_full["qval"] < 0.05, "g"].astype(str).tolist()
                elif "padj" in svg_full.columns and "g" in svg_full.columns:
                    genes = svg_full.loc[svg_full["padj"] < 0.05, "g"].astype(str).tolist()
                elif "gene" in svg_full.columns:
                    genes = svg_full["gene"].astype(str).tolist()
                elif "g" in svg_full.columns:
                    genes = svg_full["g"].astype(str).tolist()
                else:
                    genes = list(map(str, svg_full.index.tolist()))

            print("[gene_selection] svg")
        else:
            raise ValueError(
                f"Unknown gene selection mode '{mode}'. "
                "Valid modes: ctg, hvg, spapros, svg"
            )

        genes = list(genes) if genes is not None else []
        selected_sources[mode] = genes
        print(f"[gene_selection] {mode}: n={len(genes)}")

    if user_genes:
        selected_sources["user_list"] = list(dict.fromkeys(user_genes))
        print(f"[gene_selection] user_list: n={len(selected_sources['user_list'])}")

    if not selected_sources:
        print("[gene_selection] Using all overlapping genes (no selection mode / no user gene list)")
        return None

    union_genes = []
    for _, genes in selected_sources.items():
        union_genes.extend(genes)
    union_genes = list(dict.fromkeys(union_genes))

    overlap = set(ad_sc.var_names).intersection(set(ad_sp.var_names))
    union_overlap = [g for g in union_genes if g in overlap]

    print(f"[gene_selection] Union before overlap filter: n={len(union_genes)}")
    print(f"[gene_selection] Union after overlap filter:  n={len(union_overlap)}")

    if len(union_overlap) == 0:
        raise ValueError(
            "No selected genes overlap between single-cell and spatial data. "
            "Check whether your gene identifiers match ad_sc.var_names / ad_sp.var_names."
        )

    return union_overlap


def run_tangram(ad_sc: object, ad_sp: object, gene_selection_mode= None, cell_label: str = 'cell_type', ensembl_col: str = "", feature_col: str = "", device_choice: str = 'cpu', multiome: bool = False, user_genes=None):
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
    genes = select_genes(ad_sc_use, ad_sp, gene_selection_mode, cell_label, user_genes=user_genes)

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



def add_dissociated_annotations_to_obs(
    ad_ge,
    pred_key: str = "tangram_ct_pred",
    label_out: str = "cell_type_dissociated",
    prob_prefix: str = "prob_",
    prob_suffix: str = "_dissociated",
    sanitize_names: bool = True,
    add_masked_prob_columns: bool = True,
    mask_obs_key: str = "cell_type",
    masked_sep: str = "__",
):
    if pred_key not in ad_ge.obsm:
        raise KeyError(f"{pred_key!r} not found in ad_ge.obsm")

    ct_pred = ad_ge.obsm[pred_key]

    # Ensure DataFrame with proper index/column names
    if not isinstance(ct_pred, pd.DataFrame):
        ct_pred = pd.DataFrame(ct_pred, index=ad_ge.obs_names)

    ct_pred = ct_pred.copy()
    ct_pred.index = ad_ge.obs_names
    ct_pred.columns = ct_pred.columns.astype(str)

    def _safe_name(x):
        x = str(x)
        if sanitize_names:
            x = re.sub(r"[^0-9a-zA-Z_]+", "_", x).strip("_")
        return x

    # highest-score cell type per spot
    ad_ge.obs[label_out] = ct_pred.idxmax(axis=1).astype("category")

    # probability/score columns per dissociated cell type
    colmap = {}
    for ct in ct_pred.columns:
        ct_name = str(ct)
        ct_safe = _safe_name(ct_name)
        out_col = f"{prob_prefix}{ct_safe}{prob_suffix}"
        ad_ge.obs[out_col] = ct_pred[ct_name].to_numpy()
        colmap[ct_name] = out_col

    ad_ge.uns["dissociated_prob_column_map"] = colmap

    # Optional: create masked columns for plotting by existing spatial annotation
    if add_masked_prob_columns:
        if mask_obs_key not in ad_ge.obs.columns:
            raise KeyError(f"{mask_obs_key!r} not found in ad_ge.obs")

        mask_series = ad_ge.obs[mask_obs_key].astype(str)
        mask_categories = sorted(mask_series.unique())

        masked_colmap = {}
        for mask_ct in mask_categories:
            mask_ct_safe = _safe_name(mask_ct)
            mask = mask_series == mask_ct

            masked_colmap[mask_ct] = {}
            for original_ct_name, prob_col in colmap.items():
                new_col = (
                    f"{prob_col}"
                    f"{masked_sep}{_safe_name(mask_obs_key)}"
                    f"{masked_sep}{mask_ct_safe}"
                )

                vals = ad_ge.obs[prob_col].to_numpy(dtype=float).copy()
                vals[~mask.to_numpy()] = np.nan
                ad_ge.obs[new_col] = vals

                masked_colmap[mask_ct][original_ct_name] = new_col

        ad_ge.uns["dissociated_masked_prob_column_map"] = masked_colmap

    return ad_ge

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--sc_path", required=True, help="Pfad zur Single-Cell h5ad")
    parser.add_argument("--sp_path", required=True, help="Pfad zur Spatial h5ad")
    parser.add_argument("--multiome", action='store_true', help='Whether the single-cell data is from a multiome experiment (default: False)')
    parser.add_argument("--outdir", required=True, help="Pfad für Ausgabe-Verzeichnis (Tangram)")
    parser.add_argument("--gene_selection_mode", nargs="*", default=None, help="One or more gene selection modes: ctg hvg spapros svg. "
         "Also accepts comma-separated values, e.g. 'ctg,hvg'.")
    parser.add_argument("--cell_label", default="cell_type")
    parser.add_argument("--ensembl_col", default="")
    parser.add_argument("--feature_col", default="")
    parser.add_argument("--device", default="cpu", choices=["cpu", "gpu"])
    parser.add_argument("--gene_list_path", default="", help="Optional path to a user-provided gene list (txt/csv/tsv).")
    parser.add_argument("--gene_list_column", default="", help="Optional column name for csv/tsv gene lists. "
         "If empty, the first column is used."
    )
    args = parser.parse_args()

    # Daten laden
    print("[tangram_cli] Lade Single-Cell AnnData:", args.sc_path)
    ad_sc = sc.read_h5ad(args.sc_path)

    print("[tangram_cli] Lade Spatial AnnData:", args.sp_path)
    ad_sp = sc.read_h5ad(args.sp_path)

    user_genes = []
    if args.gene_list_path:
        user_genes = load_user_genes(args.gene_list_path, args.gene_list_column)

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
        user_genes=user_genes,
    )

    # Ergebnis speichern
    print("[tangram_cli] Schreibe Ergebnis nach:", args.outdir)
    print(ad_ge)
    print(ad_ge.obsm["tangram_ct_pred"])
    ad_ge.obsm["tangram_ct_pred"] = ad_ge.obsm["tangram_ct_pred"].rename(columns=str)
    ad_ge = add_dissociated_annotations_to_obs(
        ad_ge=ad_ge,
        pred_key="tangram_ct_pred",
        label_out="cell_type_dissociated",
        prob_prefix="prob_",
        prob_suffix="_dissociated",
        sanitize_names=True,
        add_masked_prob_columns=True,
        mask_obs_key="cell_type",
    )
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


