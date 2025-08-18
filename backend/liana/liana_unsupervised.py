import decoupler as dc
import pandas as pd
import scanpy as sc


def pathway_activities(
    adata: sc.AnnData,
    pathway_net: pd.DataFrame,
    return_scores: bool = False,
) -> dict | None:
    # Use multivariate linear model to estimate activity
    dc.mt.mlm(
        data=adata,
        net=pathway_net,
        verbose=True,
        raw=False,
    )

    mlm_scores = adata.obsm["score_mlm"]
    mlm_padj = adata.obsm["padj_mlm"]

    if return_scores:
        return {"score_mlm": mlm_scores, "padj_mlm": mlm_padj}
    else:
        adata.obsm["pathway_activity_score_mlm"] = mlm_scores
        adata.obsm["pathway_activity_padj_mlm"] = mlm_padj
        del adata.obsm["score_mlm"], adata.obsm["padj_mlm"]
