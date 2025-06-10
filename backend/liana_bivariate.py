import decoupler as dc
import liana as li
import mygene
import numpy as np
import pandas as pd
import scanpy as sc
from matplotlib import pyplot as plt


def convert_ensembl_to_symbol(ensembl_series, species="human"):
    mg = mygene.MyGeneInfo()
    unique_ids = ensembl_series.unique().tolist()
    results = mg.querymany(
        unique_ids, scopes="ensembl.gene", fields="symbol", species=species
    )
    id_to_symbol = {}
    for item in results:
        if "symbol" in item:
            id_to_symbol[item["query"]] = item["symbol"]
        else:
            id_to_symbol[item["query"]] = None
    return ensembl_series.map(id_to_symbol)


if __name__ == "__main__":
    # Source: https://liana-py.readthedocs.io/en/latest/notebooks/bivariate.html
    # NOTE: Make sure that liana >= 1.0.0 is installed

    # Data loading
    dataset_path = "/nfs/data3/mopitas/mapra/datasets/GSM6592049_M2/GSM6592049_M2.h5ad"
    adata = sc.read(dataset_path)
    adata = adata[adata.obs["in_tissue"] == 1]

    # Normalization
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

    # Plot spatial data and cell type annotations
    sc.pl.spatial(
        adata,
        color=[None, "cell_type"],
        size=1.3,
        palette="Set1",
        library_id="GSM6592049_M2",
    )

    plot, bw_df = li.utils.query_bandwidth(
        coordinates=adata.obsm["spatial"], start=0, end=500, interval_n=20
    )
    # plot.show()
    bw = int(np.round(bw_df[bw_df["neighbours"] >= 6].iloc[0, 0]))

    # Construct spatial connectivities
    # Bandwidth may be set based on query_bandwidth function
    li.ut.spatial_neighbors(
        adata, bandwidth=bw, cutoff=0.1, kernel="gaussian", set_diag=True
    )
    connectivity_plot = li.pl.connectivity(adata, idx=0, size=1.3, figure_size=(6, 5))
    # connectivity_plot.show()

    ensembl_ids = adata.var.index
    gene_symbols = convert_ensembl_to_symbol(ensembl_ids)

    adata.var["gene_ids"] = adata.var.index
    adata.var.index = gene_symbols
    # NOTE: required to fix common bug (MergeError)
    # https://github.com/saezlab/liana-py/issues/143
    adata.var.index.name = None

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

    sc.pl.spatial(
        lrdata,
        color=["TIMP2^ITGB1", "COL1A1^SDC1"],
        size=1.4,
        vmax=1,
        cmap="magma",
        library_id="GSM6592049_M2",
    )

    sc.pl.spatial(
        adata,
        color=["TIMP2", "ITGB1", "COL1A1", "SDC1"],
        size=1.4,
        ncols=2,
        library_id="GSM6592049_M2",
    )

    # Permutation-based p-values
    sc.pl.spatial(
        lrdata,
        layer="pvals",
        color=["TIMP2^ITGB1", "COL1A1^SDC1"],
        size=1.4,
        cmap="magma_r",
        library_id="GSM6592049_M2",
    )

    # Local Categories for ligand/receptor expression
    # +1 = high-high
    # 0 = high-low / low-high
    # -1 = low-low
    sc.pl.spatial(
        lrdata,
        layer="cats",
        color=["TIMP2^ITGB1", "COL1A1^SDC1"],
        size=1.4,
        cmap="coolwarm",
        library_id="GSM6592049_M2",
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

    lr_loadings = li.ut.get_variable_loadings(lrdata, varm_key="NMF_H").set_index(
        "index"
    )
    factor_scores = li.ut.get_factor_scores(lrdata, obsm_key="NMF_W")

    nmf = sc.AnnData(
        X=lrdata.obsm["NMF_W"],
        obs=lrdata.obs,
        var=pd.DataFrame(index=lr_loadings.columns),
        uns=lrdata.uns,
        obsm=lrdata.obsm,
    )

    sc.pl.spatial(
        nmf, color=[*nmf.var.index, None], size=1.4, ncols=2, library_id="GSM6592049_M2"
    )

    # Beyond Ligand-Receptors

    # NOTE: we don't have compositions in obsm
    # would require running cell2location, so ignore for now

    # Extract Cell type Composition
    comps = li.ut.obsm_to_adata(adata, "compositions")
    # check key cell types
    sc.pl.spatial(comps, color=["vSMCs", "CM", "Endo", "Fib"], size=1.3, ncols=2)
