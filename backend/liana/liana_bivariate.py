from itertools import product

import decoupler as dc
import liana as li
import pandas as pd
import scanpy as sc
from mudata import MuData


def ligand_receptor_relationships(
    adata: sc.AnnData, return_scores: bool = False
):
    # HACK: Setting index.name to None is required to fix
    # common bug (MergeError). If everything goes down the
    # drain because of that we might need to come back.
    # https://github.com/saezlab/liana-py/issues/143
    adata.var.index.name = None
    adata.var.index = adata.var.index.str.upper()

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
    # Make sure the GRN column names are ok
    net = net.rename(
        columns={
            old: new
            for old, new in zip(
                net.columns[:3], ["source", "target", "weight"]
            )
        }
    )

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
