import cell2location
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc

if __name__ == "__main__":
    results_folder = "../data/cell2location-results/"

    ref_run_name = f"{results_folder}/reference_signatures"
    run_name = f"{results_folder}/cell2location_map"

    adata_vis = sc.read_h5ad("/nfs/data3/mopitas/mapra/datasets/GSM6592052_M5/GSM6592052_M5.h5ad")
    adata_vis = adata_vis[adata_vis.obs["in_tissue"] == 1]
    adata_vis.obs["sample"] = list(adata_vis.uns["spatial"].keys())[0]

    adata_vis.var["MT_gene"] = [
        gene.startswith("MT-") for gene in adata_vis.var["feature_name"]
    ]

    adata_vis.obsm["MT"] = adata_vis[
        :, adata_vis.var["MT_gene"].values
    ].X.toarray()
    adata_vis = adata_vis[:, ~adata_vis.var["MT_gene"].values]

    adata_ref = sc.read("/nfs/data3/mopitas/mapra/datasets/scRNA/Wu_annotated.h5ad")

    adata_vis.var["ensembl_id"] = adata_vis.var.index
    # rename 'GeneID-2' as necessary for your data
    adata_vis.var.set_index("feature_name", drop=True, inplace=True)

    # delete unnecessary raw slot (to be removed in a future version of the tutorial)
    del adata_ref.raw

    selected = cell2location.utils.filtering.filter_genes(
        adata_ref,
        cell_count_cutoff=5,
        cell_percentage_cutoff2=0.03,
        nonz_mean_cutoff=1.12,
    )

    # filter the object
    adata_ref = adata_ref[:, selected].copy()

    # prepare anndata for the regression model
    cell2location.models.RegressionModel.setup_anndata(
        adata=adata_ref,
        # 10X reaction / sample / batch
        batch_key="sample",
        # cell type, covariate used for constructing signatures
        labels_key="cell_subclass",
        # multiplicative technical effects (platform, 3' vs 5', donor effect)
        # categorical_covariate_keys=[""],
    )

    # create the regression model
    mod = cell2location.models.RegressionModel(adata_ref)

    mod.train(max_epochs=250)

    mod.plot_history(20)
    # plt.show()

    adata_ref = mod.export_posterior(
        adata_ref,
        sample_kwargs={
            "num_samples": 1000,
            "batch_size": 2500,
        },
    )

    if "means_per_cluster_mu_fg" in adata_ref.varm.keys():
        inf_aver = adata_ref.varm["means_per_cluster_mu_fg"][
            [
                f"means_per_cluster_mu_fg_{i}"
                for i in adata_ref.uns["mod"]["factor_names"]
            ]
        ].copy()
    else:
        inf_aver = adata_ref.var[
            [
                f"means_per_cluster_mu_fg_{i}"
                for i in adata_ref.uns["mod"]["factor_names"]
            ]
        ].copy()
    inf_aver.columns = adata_ref.uns["mod"]["factor_names"]
    inf_aver.iloc[0:5, 0:5]

    # Cell2location spatial mapping

    intersect = np.intersect1d(adata_vis.var_names, inf_aver.index)
    adata_vis = adata_vis[:, intersect].copy()
    inf_aver = inf_aver.loc[intersect, :].copy()

    # prepare anndata for cell2location model
    cell2location.models.Cell2location.setup_anndata(
        adata=adata_vis, batch_key="sample"
    )

    mod = cell2location.models.Cell2location(
        adata_vis,
        cell_state_df=inf_aver,
        # the expected average cell abundance: tissue-dependent
        # hyper-prior which can be estimated from paired histology:
        N_cells_per_location=30,
        # hyperparameter controlling normalisation of
        # within-experiment variation in RNA detection:
        detection_alpha=20,
    )

    mod.train(
        max_epochs=30000,
        # train using full data (batch_size=None)
        batch_size=None,
        # use all data points in training because
        # we need to estimate cell abundance at all locations
        train_size=1,
    )

    # plot ELBO loss history during training, removing first 100 epochs from the plot
    mod.plot_history(1000)
    plt.legend(labels=["full data training"])
    # plt.show()

    adata_vis = mod.export_posterior(
        adata_vis,
        sample_kwargs={
            "num_samples": 1000,
            "batch_size": mod.adata.n_obs,
        },
    )

    adata_vis.write("/nfs/data3/mopitas/mapra/user/noah/cell2location_adata.h5ad")
