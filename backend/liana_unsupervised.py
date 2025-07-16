import decoupler as dc
import liana as li
import pandas as pd
import scanpy as sc

# Helper functions to create MISTy objects
from liana.method import MistyData, genericMistyData, lrMistyData

# Pre-defined single view models
from liana.method.sp import LinearModel, RandomForestModel, RobustLinearModel


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


if __name__ == "__main__":
    # Dataset loading
    dataset_path = "/home/noah/Downloads/tangram_adata.h5ad"
    adata = sc.read_h5ad(dataset_path)
    adata = adata[adata.obs["in_tissue"] == 1]
    sc.pp.filter_cells(adata, min_counts=1)

    # Normalization
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    adata.var.index.name = None
    adata.var.index = adata.var.index.str.upper()

    # comps = li.ut.obsm_to_adata(adata, "tangram_ct_pred")
    # sc.pl.spatial(
    #     comps,
    #     color=[
    #         "T-cells",
    #         "PVL",
    #         "Endothelial",
    #         "B-cells",
    #     ],
    #     size=1.3,
    #     ncols=2,
    #     alpha_img=0,
    #     library_id="GSM6592052_M5",
    # )

    pathway_net = dc.op.progeny(organism="human", top=500)
    pathway_activities(adata, pathway_net)

    # sc.pl.spatial(
    #     adata,
    #     color=["p53_pathway_activity_score"],
    #     cmap="RdBu_r",
    #     size=1.3,
    #     library_id="GSM6592052_M5",
    # )
    #
    # sc.pl.spatial(
    #     adata,
    #     color=["p53_pathway_activity_padj"],
    #     cmap="magma_r",
    #     size=1.3,
    #     library_id="GSM6592052_M5",
    # )

    # # Formatting and running MISTy
    # # MISTy gives high level summary statistics (aka global relationships)
    #
    # misty = genericMistyData(
    #     intra=comps,
    #     extra=acts_progeny,
    #     cutoff=0.05,
    #     bandwidth=100,
    #     n_neighs=6,
    #     add_juxta=False,
    #     add_para=False,
    # )
    #
    # misty(model=RandomForestModel, n_jobs=-1, verbose=True)
    #
    # misty(model=LinearModel, k_cv=10, seed=1337, verbose=True)
    #
    # # Build Custom Misty Views
    #
    # # get TF prior knowledge
    # net = dc.op.collectri(organism="human")
    #
    # # Estimate activities
    # dc.mt.ulm(data=adata, net=net, verbose=True, raw=False)
    #
    # # extract activities
    # acts_tfs = li.ut.obsm_to_adata(adata, "score_ulm")
    #
    # li.ut.spatial_neighbors(
    #     acts_tfs, cutoff=0.1, bandwidth=200, set_diag=False
    # )
    #
    # # li.pl.connectivity(acts_tfs, idx=0, figure_size=(6, 5))
