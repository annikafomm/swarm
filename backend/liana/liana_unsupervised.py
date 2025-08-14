import decoupler as dc
import pandas as pd
import scanpy as sc


def pathway_activities(
    adata: sc.AnnData,
    pathway_net: pd.DataFrame,
    return_scores: bool = False,
):
    # Use multivariate linear model to estimate activity
    dc.mt.mlm(
        data=adata,
        net=pathway_net,
        verbose=True,
        raw=False,
    )

    mlm_scores = adata.obsm["score_mlm"].add_suffix("_pathway_activity_score")
    mlm_padj = adata.obsm["padj_mlm"].add_suffix("_pathway_activity_padj")
    scores = mlm_scores.join(mlm_padj)

    if return_scores:
        return scores
    else:
        adata.obs = adata.obs.join(scores)
