from itertools import product

import decoupler as dc
import liana as li
import pandas as pd
import scanpy as sc
from mudata import MuData
import matplotlib
matplotlib.use('Agg')

CAT_REMAP = {-1: "high-low / low-high", 0: "low-low", 1: "high-high"}


def _index_by_first_column(df: pd.DataFrame) -> pd.DataFrame:
    """Re-index a liana loadings/scores frame on its identifier column.

    These frames come back from `li.ut.get_variable_loadings` / `li.ut.get_factor_scores`
    already `reset_index()`-ed, with the identifier in the first column and `Factor1..N`
    after it. The two helpers disagree about that column's name:
    `get_variable_loadings` clears `df.index.name` before resetting, so it always yields
    `"index"`, while `get_factor_scores` does not — it yields whatever `adata.obs.index.name`
    happens to be. Hardcoding `.set_index("index")` therefore worked only for datasets whose
    obs index is unnamed, and raised `KeyError: "None of ['index'] are in the columns"` for
    any dataset that names it (e.g. the heart multiome builtin, whose index is `spot_id`).

    Selecting the first column instead is correct for both, since the factor columns are
    always appended after the identifier.
    """
    if df.columns.empty:
        raise ValueError("expected an identifier column plus factor columns, got none")
    return df.set_index(df.columns[0])


def ligand_receptor_relationships(
    adata: sc.AnnData, return_scores: bool = False, resource_name: str = "consensus"
) -> dict | None:
    # Bivariate Ligand-Receptor Relationships
    # Parameters from tutorial
    #
    # `resource_name` must match the organism's symbol convention: "consensus" is
    # human-symbol based, "mouseconsensus" is mouse. Passing the wrong one silently yields
    # (almost) no matching LR pairs rather than an error. run_liana picks it from the
    # organism, so symbols never need to be rewritten to fit the resource — see the note in
    # calc_liana.run_liana about the upper-casing this replaced.
    # `li.rs.show_resources()` lists the rest; other organisms go via
    # `li.rs.get_hcop_orthologs()` + `li.rs.translate_resource()`.
    lrdata = li.mt.bivariate(
        adata,
        resource_name=resource_name,
        local_name="cosine",  # Name of the function
        global_name="morans",  # Name global function
        n_perms=100,  # Number of permutations to calculate a p-value
        mask_negatives=False,  # Whether to mask LowLow/NegativeNegative interactions
        add_categories=True,  # Whether to add local categories to the results
        nz_prop=0.2,  # Minimum expr. proportion for ligands/receptors and their subunits
        use_raw=False,
        verbose=True,
    )
    # Extract ligand-receptor scores
    lr_relationships = lrdata.var.index.to_list()
    lr_cos_sim_mat = lrdata.X.toarray()
    lr_pval_mat = lrdata.layers["pvals"].toarray()
    lr_cat_mat = lrdata.layers["cats"].toarray()
    cos_sim_df = pd.DataFrame(
        lr_cos_sim_mat,
        columns=lr_relationships,
        index=adata.obs.index,
    )
    pval_df = pd.DataFrame(
        lr_pval_mat,
        columns=lr_relationships,
        index=adata.obs.index,
    )
    cat_df = pd.DataFrame(
        lr_cat_mat,
        columns=lr_relationships,
        index=adata.obs.index,
    )
    cat_df = cat_df.replace(CAT_REMAP)

    global_score_columns = ["mean", "std", "morans"]
    global_scores = lrdata.var[global_score_columns]
    global_scores = global_scores.rename(
        columns={
            col: (
                f"cosine_similarity_{col}"
                if col in ["mean", "std"]
                else f"ligand_receptor_{col}"
            )
            for col in global_score_columns
        }
    )

    # Identify intercellular patterns
    li.multi.nmf(
        lrdata,
        n_components=None,
        inplace=True,
        random_state=0,
        max_iter=200,
        verbose=True,
    )

    # Loadings have shape ligand receptor pairs x number of factors.
    # Don't think we can use them in a straightforward way, but
    # could come up with something later.
    lr_loadings = _index_by_first_column(
        li.ut.get_variable_loadings(lrdata, varm_key="NMF_H")
    )

    factor_scores = _index_by_first_column(
        li.ut.get_factor_scores(lrdata, obsm_key="NMF_W")
    )

    if return_scores:
        local_scores = {
            "cosine_similarity": cos_sim_df,
            "p_value": pval_df,
            "category": cat_df,
        }
        return {
            "local_scores": local_scores,
            "global_scores": global_scores,
            "NMF_factors": factor_scores,
        }
    else:
        # HACK: It would make more sense to save the dfs to obsm, but the
        # number of columns we can write in a df in obsm (or uns) apparently is
        # limited. If exceeded, write_h5ad fails with `OSError: Unable to
        # synchronously create attribute (object header message is too large)`
        # https://github.com/h5py/h5py/issues/1053
        #
        # As a workaround we can reconstruct the df after reading in the adata
        # from h5ad like so:
        # adata.obsm["ligand_receptor_cosine_similarity"] = pd.DataFrame(
        #     adata.obsm["ligand_receptor_cosine_similarity"],
        #     columns=adata.uns["liana_columns"]["ligand_receptor"],
        #     index=adata.obs_names
        # )
        adata.obsm["ligand_receptor_cosine_similarity"] = cos_sim_df.to_numpy()
        adata.obsm["ligand_receptor_p_value"] = pval_df.to_numpy()
        adata.obsm["ligand_receptor_category"] = cat_df.to_numpy()
        adata.uns["liana_columns"][
            "ligand_receptor"
        ] = cos_sim_df.columns.tolist()
        adata.obsm["ligand_receptor_NMF_factors"] = factor_scores
        adata.uns["ligand_receptor_global_scores"] = global_scores


def cell_comp_tf_activity_similarity(
    adata: sc.AnnData,
    net: pd.DataFrame,
    cell_comp_obsm_key: str = "tangram_ct_pred",
    return_scores: bool = False,
) -> dict | None:
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

    # Cosine similarities
    comp_tf_interactions = bdata.var.index.to_list()
    comp_tf_cos_sim_mat = bdata.X.toarray()
    cos_sim_df = pd.DataFrame(
        comp_tf_cos_sim_mat,
        columns=comp_tf_interactions,
        index=adata.obs.index,
    )

    comp_tf_cat_mat = bdata.layers["cats"].toarray()
    cat_df = pd.DataFrame(
        comp_tf_cat_mat,
        columns=comp_tf_interactions,
        index=adata.obs.index,
    )
    cat_df = cat_df.replace(CAT_REMAP)

    global_score_columns = ["mean", "std"]
    global_scores = bdata.var[global_score_columns]
    global_scores = global_scores.rename(
        columns={
            col: f"cosine_similarity_{col}" for col in global_score_columns
        }
    )

    # ULM scores for TFs
    ulm_scores = adata.obsm["score_ulm"]
    ulm_padj = adata.obsm["padj_ulm"]
    if return_scores:
        return {
            "cosine_similarity": cos_sim_df,
            "category": cat_df,
            "global_scores": global_scores,
            "score_ulm": ulm_scores,
            "padj_ulm": ulm_padj,
        }
    else:
        adata.obsm["cell_comp_tf_activity_cosine_similarity"] = (
            cos_sim_df.to_numpy()
        )
        adata.obsm["cell_comp_tf_activity_category"] = cat_df.to_numpy()
        adata.uns["liana_columns"][
            "cell_comp_tf_activity"
        ] = cos_sim_df.columns.tolist()
        adata.obsm["tf_activity_score_ulm"] = ulm_scores
        adata.obsm["tf_activity_padj_ulm"] = ulm_padj
        del adata.obsm["score_ulm"], adata.obsm["padj_ulm"]
        adata.uns["cell_comp_tf_activity_global_scores"] = global_scores
