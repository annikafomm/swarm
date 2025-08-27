# like Noah's run_liana.py file, but slight structure changes

import decoupler as dc
import pandas as pd
import scanpy as sc
from liana_folder.liana_bivariate import (
    cell_comp_tf_activity_similarity,
    ligand_receptor_relationships,
)
from liana_folder.liana_unsupervised import pathway_activities

def run_liana(adata, grn=None, pathway_net=None, cell_comp_key="tangram_ct_pred"):
    adata.uns["liana_columns"] = {}

    if grn is None:
        grn = dc.op.collectri()
    else:
        grn = pd.read_csv(grn)

    if pathway_net is None:
        pathway_net = dc.op.progeny()
    else:
        pathway_net = pd.read_csv(pathway_net)

    adata.var.index.name = None
    adata.var.index = adata.var.index.str.upper()

    ligand_receptor_relationships(adata)
    cell_comp_tf_activity_similarity(adata, grn, cell_comp_obsm_key=cell_comp_key)
    pathway_activities(adata, pathway_net)

    return adata
