"""
Source: https://liana-py.readthedocs.io/en/latest/notebooks/bivariate.html
NOTE: Make sure that liana >= 1.0.0 is installed
"""

import argparse

import decoupler as dc
import pandas as pd
import scanpy as sc
import squidpy as sq
from liana_bivariate import (
    cell_comp_tf_activity_similarity,
    ligand_receptor_relationships,
)
from liana_unsupervised import pathway_activities

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--adata", type=str, required=True)
    parser.add_argument("--grn", type=str, default=None)
    parser.add_argument("--pathway-net", type=str, default=None)
    parser.add_argument("--output", type=str, required=True)

    args = parser.parse_args()

    adata = sc.read_h5ad(args.adata)

    if args.grn is not None:
        grn = pd.read_csv(args.grn)
    else:
        grn = dc.op.collectri()

    if args.pathway_net is not None:
        pathway_net = pd.read_csv(args.pathway_net)
    else:
        pathway_net = dc.op.progeny()

    ligand_receptor_relationships(adata)

    cell_comp_tf_activity_similarity(
        adata,
        grn,
    )

    pathway_activities(adata, pathway_net)

    adata.write_h5ad(args.output)
