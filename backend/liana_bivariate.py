from itertools import product

import decoupler as dc
import liana as li
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt
from mudata import MuData


def ligand_receptor_relationships(
    adata: sc.AnnData, return_scores: bool = False
):
    # Bivariate Ligand-Receptor Relationships
    # Parameters from tutorial
    lrdata = li.mt.bivariate(
        adata,
        resource_name="consensus",  # NOTE: uses HUMAN gene symbols!
        local_name="cosine",  # Name of the function
        global_name="morans",  # Name global function
        n_perms=100,  # Number of permutations to calculate a p-value
        mask_negatives=False,  # Whether to mask LowLow/NegativeNegative interactions
        add_categories=True,  # Whether to add local categories to the results
        nz_prop=0.2,  # Minimum expr. proportion for ligands/receptors and their subunits
        use_raw=False,
        verbose=True,
    )
    # Add ligand-receptor scores to adata.obs
    lr_relationships = lrdata.var.index
    lr_cos_sim_mat = lrdata.X.T.toarray()
    lr_pval_mat = lrdata.layers["pvals"].T.toarray()
    lr_cat_mat = lrdata.layers["cats"].T.toarray()
    local_scores = pd.DataFrame()
    for i, lr_relationship in enumerate(lr_relationships):
        local_scores[f"{lr_relationship}_cosine-similarity"] = lr_cos_sim_mat[
            i
        ]
        local_scores[f"{lr_relationship}_p-value"] = lr_pval_mat[i]
        local_scores[f"{lr_relationship}_category"] = lr_cat_mat[i]

    global_scores = lrdata.var[["mean", "std", "morans"]]

    # Identify intercellular patterns
    li.multi.nmf(
        lrdata,
        n_components=None,
        inplace=True,
        random_state=0,
        max_iter=200,
        verbose=True,
    )

    # Loadings have shape ligand receptor pairs x number of factors
    # Don't think we can use them
    lr_loadings = li.ut.get_variable_loadings(
        lrdata, varm_key="NMF_H"
    ).set_index("index")

    factor_scores = li.ut.get_factor_scores(
        lrdata, obsm_key="NMF_W"
    ).set_index("index")

    factor_scores = factor_scores.rename(
        columns={
            f"Factor{i}": f"NMF_Factor_{i}"
            for i, _ in enumerate(factor_scores.columns, 1)
        }
    )

    # adata.obs = adata.obs.join(factor_scores)
    if return_scores:
        return (local_scores, global_scores, factor_scores)
    else:
        adata.obs = (
            adata.obs.join(local_scores)
            .join(global_scores)
            .join(factor_scores)
        )


def cell_comp_tf_activity_similarity(
    adata: sc.AnnData,
    net: pd.DataFrame,
    cell_comp_obsm_key: str = "tangram_ct_pred",
    return_scores: bool = False,
):
    # Extract Cell type Composition
    comps = li.ut.obsm_to_adata(adata, cell_comp_obsm_key)

    dc.mt.ulm(adata, net=net, raw=False, verbose=True)

    est = li.ut.obsm_to_adata(adata, "score_ulm")
    est.var["cv"] = est.X.std(axis=0) / est.X.mean(axis=0)
    top_tfs = (
        est.var.sort_values("cv", ascending=False, key=abs).head(50).index
    )

    mdata = MuData({"tf": est, "comps": comps})
    mdata.obsp = adata.obsp
    mdata.uns = adata.uns
    mdata.obsm = adata.obsm

    interactions = list(product(comps.var.index, top_tfs))

    bdata = li.mt.bivariate(
        mdata,
        x_mod="comps",
        y_mod="tf",
        x_transform=sc.pp.scale,
        y_transform=sc.pp.scale,
        local_name="cosine",
        interactions=interactions,
        mask_negatives=True,
        add_categories=True,
        x_use_raw=False,
        y_use_raw=False,
        xy_sep="<->",
        x_name="celltype",
        y_name="tf",
    )

    scores = pd.DataFrame()
    # Cosine similarities
    comp_tf_interactions = bdata.var.index
    comp_tf_cos_sim_mat = bdata.X.T.toarray()
    for interaction, cos_sim_vec in zip(
        comp_tf_interactions, comp_tf_cos_sim_mat
    ):
        scores[interaction] = cos_sim_vec

    # ULM scores for TFs
    tfs = mdata.mod["tf"].var.index
    ulm_score_mat = mdata.mod["tf"].X.T
    for tf, ulm_score_vec in zip(tfs, ulm_score_mat):
        scores[tf] = ulm_score_vec

    if return_scores:
        return scores
    else:
        adata.obs = adata.obs.join(scores)


if __name__ == "__main__":
    # Source: https://liana-py.readthedocs.io/en/latest/notebooks/bivariate.html
    # NOTE: Make sure that liana >= 1.0.0 is installed

    # Data loading
    # dataset_path = "/nfs/data3/mopitas/mapra/datasets/GSM6592049_M2/GSM6592049_M2.h5ad"
    dataset_path = "/home/noah/Downloads/tangram_adata.h5ad"
    c2l_path = "/home/noah/Downloads/cell2location_adata.h5ad"
    adata = sc.read(dataset_path)
    c2l_adata = sc.read(c2l_path)
    adata.obsm["c2l_compositions"] = c2l_adata.obsm["q05_cell_abundance_w_sf"]

    # Only look at spots in tissue
    adata = adata[adata.obs["in_tissue"] == 1]
    # Make sure that all cells have counts, else decoupler's ULM throws them
    # out with just a warning, causing annoying shape mismatches
    sc.pp.filter_cells(adata, min_counts=1)

    # Normalization from tutorial
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Spatial connectivity graph
    plot, bw_df = li.utils.query_bandwidth(
        coordinates=adata.obsm["spatial"], start=0, end=500, interval_n=20
    )
    bw = int(np.round(bw_df[bw_df["neighbours"] >= 6].iloc[0, 0]))

    # Bandwidth may be set based on query_bandwidth function
    li.ut.spatial_neighbors(
        adata, bandwidth=bw, cutoff=0.1, kernel="gaussian", set_diag=True
    )

    # NOTE: Setting index.name to None is required to fix
    # common bug (MergeError)
    # https://github.com/saezlab/liana-py/issues/143
    adata.var.index.name = None
    adata.var.index = adata.var.index.str.upper()

    ligand_receptor_relationships(adata)

    # net = dc.op.collectri()
    net = pd.read_csv(
        "/home/noah/Downloads/genie3_BRCA_mrn.top_100k.csv",
    )
    net = net.rename(
        columns={
            old: new
            for old, new in zip(net.columns, ["source", "target", "weight"])
        }
    )

    cell_comp_tf_activity_similarity(adata, net, "c2l_compositions")

    adata.write("/home/noah/Downloads/liana_cell2location_adata.h5ad")
