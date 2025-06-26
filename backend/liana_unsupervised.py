import decoupler as dc
import liana as li
import scanpy as sc

# Helper functions to create MISTy objects
from liana.method import MistyData, genericMistyData, lrMistyData

# Pre-defined single view models
from liana.method.sp import LinearModel, RandomForestModel, RobustLinearModel

if __name__ == "__main__":
    # Dataset loading
    dataset_path = "/home/noah/Downloads/cell2location_adata.h5ad"
    adata = sc.read_h5ad(dataset_path)
    adata = adata[adata.obs["in_tissue"] == 1]

    # Normalization
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    comps = li.ut.obsm_to_adata(adata, "means_cell_abundance_w_sf")

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

    # NOTE: is a GRN
    progeny = dc.op.progeny(organism="human", top=500)

    # use multivariate linear model to estimate activity
    adata = dc.mt.mlm(
        data=adata,
        net=progeny,
        verbose=True,
        raw=False,
    )

    acts_progeny = li.ut.obsm_to_adata(adata, "score_mlm")

    # sc.pl.spatial(
    #     acts_progeny,
    #     color=["PI3K"],
    #     cmap="RdBu_r",
    #     size=1.3,
    #     library_id="GSM6592052_M5",
    # )

    # Formatting and running MISTy

    misty = genericMistyData(
        intra=comps,
        extra=acts_progeny,
        cutoff=0.05,
        bandwidth=100,
        n_neighs=6,
        add_juxta=False,
        add_para=False,
    )

    misty(model=RandomForestModel, n_jobs=-1, verbose=True)

    misty(model=LinearModel, k_cv=10, seed=1337, verbose=True)

    # Build Custom Misty Views

    # get TF prior knowledge
    # NOTE: is a GRN
    net = dc.op.collectri(organism="human")

    # Estimate activities
    dc.mt.ulm(data=adata, net=net, verbose=True, raw=False)

    # extract activities
    acts_tfs = li.ut.obsm_to_adata(adata, "score_ulm")

    li.ut.spatial_neighbors(
        acts_tfs, cutoff=0.1, bandwidth=200, set_diag=False
    )

    # li.pl.connectivity(acts_tfs, idx=0, figure_size=(6, 5))
