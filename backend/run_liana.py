import argparse

import decoupler as dc
import pandas as pd
import scanpy as sc
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

    args = parser.parse_args()

    adata = sc.read_h5ad(args.adata)

    local_lr_scores, global_lr_scores, factor_scores = (
        ligand_receptor_relationships(adata, return_scores=True)
    )

    grn = (
        pd.read_csv(args.grn) if args.grn is not None else dc.op.collectri(),
    )
    comp_tf_scores = cell_comp_tf_activity_similarity(
        adata,
        grn,
        return_scores=True,
    )

    pathway_net = (
        pd.read_csv(args.pathway_net)
        if args.pathway_net is not None
        else dc.op.progeny()
    )
    pathway_scores = pathway_activities(adata, pathway_net, return_scores=True)
