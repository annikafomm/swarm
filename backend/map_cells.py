#!/usr/bin/env python3
# coding: utf-8

import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import tangram as tg
import torch
import gene_selection

def load_data(sc_path, sp_path):
    ad_sp = sc.read_h5ad(sp_path)
    ad_sc = sc.read_h5ad(sc_path)
    return ad_sc, ad_sp

def select_genes(ad_sc, ad_sp, selection_mode):
    if selection_mode == "ctg":
        genes = gene_selection.ctg(ad_sc, "cell_subclass")
        print("ctg")
    elif selection_mode == "hvg":
        genes = gene_selection.hvg(ad_sc)
        print("hvg")
    elif selection_mode == "spapros":
        sc.pp.highly_variable_genes(ad_sc, flavor="seurat", n_top_genes=2000)
        genes = gene_selection.spapros(ad_sc, "cell_subclass")
        print("spapros")
    elif selection_mode == "svg":
        ad_sp.raw = ad_sp.copy()
        genes = gene_selection.svg(ad_sp)
        print("svg")
    else:
        print("No gene selection performed. All overlapping genes will be used.")
        genes = None 
        print("none")
    return genes

def main(args):
    ad_sc, ad_sp = load_data(args.sc_path, args.sp_path)

    # Gene selection
    genes = select_genes(ad_sc, ad_sp, args.gene_selection)

    # Preprocessing
    tg.pp_adatas(ad_sc, ad_sp, genes=genes)

    # Mapping
    ad_map = tg.map_cells_to_space(ad_sc, ad_sp, device="cuda:0" if torch.cuda.is_available() else "cpu")

    # Project gene expression
    ad_ge = tg.project_genes(ad_map, ad_sc)
    ad_ge.obsm["spatial"] = ad_sp.obsm["spatial"]

    # Project cell annotations
    tg.project_cell_annotations(ad_map, ad_ge, annotation="cell_subclass")
    ad_ge.write(args.output_dir + "/ad_cell_type.h5ad")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Tangram cell mapping pipeline.")
    parser.add_argument('--sc_path', type=str, required=True, help="Path to the single-cell h5ad file.")
    parser.add_argument('--sp_path', type=str, required=True, help="Path to the spatial h5ad file.")
    parser.add_argument('--output_dir', type=str, required=True, help="Path to the output directory.")
    parser.add_argument('--gene_selection', type=str, choices=['ctg', 'hvg', 'spapros', 'svg'],
                        default=None, help="Optional gene selection: ctg, hvg, spapros, svg. Default: all genes.")

    args = parser.parse_args()
    main(args)
