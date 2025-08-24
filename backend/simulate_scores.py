"""
This script exists for debugging purposes.

Simulate all scores based on an anndata.
The output is an anndata with entries for all scores,
but the numbers will be randomly sampled from uniform
distributions whose bounds are based on real data.
"""

import argparse
import os

import numpy as np
import pandas as pd
import scanpy as sc

np.random.seed(777)


def sim_genewise_scores(adata: sc.AnnData):
    """
    - .uns["moranI"]
    - .uns["gearyC"]
    """
    num_genes = adata.shape[1]
    adata.uns["gearyC"] = pd.DataFrame(
        {
            "C": np.random.uniform(0.4, 1.7, size=num_genes),
            "pval_norm": np.random.uniform(0.0, 0.5, size=num_genes),
            "var_norm": np.random.uniform(0.0, 2e-4, size=num_genes),
            "pval_norm_fdr_bh": np.random.uniform(0.0, 0.5, size=num_genes),
        },
        index=adata.var.index,
    )
    adata.uns["moranI"] = pd.DataFrame(
        {
            "I": np.random.uniform(-0.03, 0.6, size=num_genes),
            "pval_norm": np.random.uniform(0.0, 0.5, size=num_genes),
            "var_norm": np.random.uniform(0.0, 2e-4, size=num_genes),
            "pval_norm_fdr_bh": np.random.uniform(0.0, 0.5, size=num_genes),
        },
        index=adata.var.index,
    )


def sim_leiden_scores(adata: sc.AnnData):
    """
    - .obs["leiden"]
    - .uns["leiden_centrality_scores"]
    - .uns["leiden_co_occurrence"]
    - .uns["leiden_nhood_enrichment"]
    """
    num_genes = adata.shape[1]
    num_clusters = 10
    num_bins = 49

    adata.obs["leiden"] = np.random.randint(0, num_clusters, adata.n_obs)
    adata.uns["leiden_centrality_scores"] = pd.DataFrame(
        {
            "degree_centrality": np.random.uniform(0.01, 0.2, size=num_genes),
            "average_clustering": np.random.uniform(0.4, 0.43, size=num_genes),
            "closeness_centrality": np.random.uniform(
                0.09, 0.32, size=num_genes
            ),
        }
    )

    interval = np.linspace(300, 15000, num_bins + 1, dtype=np.float32)
    occ = np.zeros((num_clusters, num_clusters, num_bins), dtype=np.float32)
    for i in range(num_clusters):
        for j in range(num_clusters):
            if i == j:
                # Strong diagonal values, log-normal distribution
                base_curve = np.random.lognormal(
                    mean=1.2, sigma=0.5, size=num_bins
                )
            else:
                # Weaker off-diagonal, sometimes zeros
                base_curve = np.random.lognormal(
                    mean=-2.0, sigma=0.3, size=num_bins
                )
                if np.random.rand() < 0.2:  # 20% chance of all zeros
                    base_curve[:] = 0.0
            # Smooth with cumulative averaging
            base_curve = np.cumsum(base_curve) / np.arange(1, num_bins + 1)
            occ[i, j, :] = base_curve.astype(np.float32)

    adata.uns["leiden_co_occurrence"] = {"occ": occ, "interval": interval}

    count = np.zeros((num_clusters, num_clusters), dtype=np.uint32)
    for i in range(num_clusters):
        for j in range(num_clusters):
            if i == j:
                # big diagonal counts
                count[i, j] = np.random.randint(500, 4000)
            else:
                # smaller off-diagonal
                count[i, j] = np.random.randint(0, 200)

    # symmetricize (like Scanpy does)
    count = (count + count.T) // 2

    # z-score: standardized version of counts
    mean = count.mean()
    std = count.std() + 1e-6
    zscore = (count - mean) / std

    adata.uns["leiden_nhood_enrichment"] = {
        "count": count,
        "zscore": zscore.astype(np.float32),
    }


def sim_genesets(adata: sc.AnnData, sponge: bool = False, num_sets: int = 5):
    genes = (
        adata.var.index.to_numpy()
        if not sponge
        else adata.var["ensemble_id"].to_numpy()
    )

    keys = np.random.choice(genes, size=num_sets, replace=False)

    genesets = {}
    for key in keys:
        size = np.random.randint(10, 50 + 1)  # random size
        genesets[key] = np.random.choice(genes, size=size, replace=False)

    adata.uns["sponge_genesets" if sponge else "genie_genesets"] = genesets


def sim_Rscores(adata: sc.AnnData):
    num_sets = 5
    num_cells = adata.shape[0]
    sim_genesets(adata, sponge=False, num_sets=num_sets)
    sim_genesets(adata, sponge=True, num_sets=num_sets)
    genie_keys = adata.uns["genie_genesets"].keys()
    sponge_keys = adata.uns["sponge_genesets"].keys()

    adata.obsm["aucell_scores_genie3"] = pd.DataFrame(
        np.random.normal(loc=0.5, scale=0.05, size=(num_cells, num_sets)),
        index=adata.obs_names,
        columns=genie_keys,
    )
    adata.obsm["aucell_scores_sponge"] = pd.DataFrame(
        np.random.normal(loc=0.5, scale=0.05, size=(num_cells, num_sets)),
        index=adata.obs_names,
        columns=sponge_keys,
    )

    adata.obsm["spongeffects_GSVA_scores_genie3"] = pd.DataFrame(
        np.random.normal(loc=-0.2, scale=0.15, size=(num_cells, num_sets)),
        index=adata.obs_names,
        columns=genie_keys,
    )
    adata.obsm["spongeffects_GSVA_scores_sponge"] = pd.DataFrame(
        np.random.normal(loc=-0.2, scale=0.15, size=(num_cells, num_sets)),
        index=adata.obs_names,
        columns=sponge_keys,
    )

    adata.obsm["spongeffects_ssGSEA_scores_genie3"] = pd.DataFrame(
        np.random.normal(loc=-0.05, scale=0.1, size=(num_cells, num_sets)),
        index=adata.obs_names,
        columns=genie_keys,
    )
    adata.obsm["spongeffects_ssGSEA_scores_sponge"] = pd.DataFrame(
        np.random.normal(loc=-0.05, scale=0.1, size=(num_cells, num_sets)),
        index=adata.obs_names,
        columns=sponge_keys,
    )

    adata.obsm["viper_scores_genie3"] = pd.DataFrame(
        np.random.normal(loc=0.0, scale=1.0, size=(num_cells, num_sets)),
        index=adata.obs_names,
        columns=genie_keys,
    )


def sim_liana_scores(adata: sc.AnnData):
    adata.uns["liana_columns"] = {}

    num_ligands = 10
    num_receptors = 8
    num_factors = 5
    lr_pairs = [
        f"lig{lig}^rec{rec}"
        for lig in range(num_ligands)
        for rec in range(num_receptors)
    ]
    adata.uns["liana_columns"]["ligand_receptor"] = np.array(lr_pairs)
    adata.obsm["ligand_receptor_cosine_similarity"] = np.random.uniform(
        0.0, 1.0, size=(adata.n_obs, num_ligands * num_receptors)
    )
    adata.obsm["ligand_receptor_p_value"] = np.random.uniform(
        0.0, 1.0, size=(adata.n_obs, num_ligands * num_receptors)
    )
    adata.obsm["ligand_receptor_category"] = np.random.randint(
        -1, 2, size=(adata.n_obs, num_ligands * num_receptors)
    )
    adata.obsm["ligand_receptor_NMF_factors"] = pd.DataFrame(
        np.random.uniform(0.0, 1.0, size=(adata.n_obs, num_factors)),
        columns=[f"Factor{f}" for f in range(1, num_factors + 1)],
        index=adata.obs.index,
    )
    adata.uns["ligand_receptor_global_scores"] = pd.DataFrame(
        {
            "cosine_similarity_mean": adata.obsm[
                "ligand_receptor_cosine_similarity"
            ].mean(axis=0),
            "cosine_similarity_std": adata.obsm[
                "ligand_receptor_cosine_similarity"
            ].std(axis=0),
            "ligand_receptor_morans": np.random.uniform(
                -0.7, 0.7, size=num_ligands * num_receptors
            ),
        },
        index=lr_pairs,
    )

    num_celltypes = 5
    num_tfs = 5
    lr_pairs = [
        f"ct{ct}<->tf{tf}"
        for ct in range(num_celltypes)
        for tf in range(num_tfs)
    ]
    adata.uns["liana_columns"]["cell_comp_tf_activity"] = np.array(lr_pairs)
    adata.obsm["cell_comp_tf_activity_cosine_similarity"] = np.random.uniform(
        0.0, 1.0, size=(adata.n_obs, num_celltypes * num_tfs)
    )
    adata.obsm["cell_comp_tf_activity_category"] = np.random.randint(
        -1, 2, size=(adata.n_obs, num_celltypes * num_tfs)
    )
    adata.uns["cell_comp_tf_activity_global_scores"] = pd.DataFrame(
        {
            "cosine_similarity_mean": adata.obsm[
                "cell_comp_tf_activity_cosine_similarity"
            ].mean(axis=0),
            "cosine_similarity_std": adata.obsm[
                "cell_comp_tf_activity_cosine_similarity"
            ].std(axis=0),
        },
        index=lr_pairs,
    )

    adata.obsm["tf_activity_score_ulm"] = pd.DataFrame(
        {
            f"tf{pw}": np.random.uniform(-10, 70, size=adata.n_obs)
            for pw in range(num_tfs)
        },
        index=adata.obs.index,
    )
    adata.obsm["tf_activity_padj_ulm"] = pd.DataFrame(
        {
            f"tf{pw}": np.random.uniform(0, 1, size=adata.n_obs)
            for pw in range(num_tfs)
        },
        index=adata.obs.index,
    )

    num_pathways = 5
    adata.obsm["pathway_activity_score_mlm"] = pd.DataFrame(
        {
            f"pathway{pw}": np.random.uniform(-10, 20, size=adata.n_obs)
            for pw in range(num_pathways)
        },
        index=adata.obs.index,
    )
    adata.obsm["pathway_activity_padj_mlm"] = pd.DataFrame(
        {
            f"pathway{pw}": np.random.uniform(0, 1, size=adata.n_obs)
            for pw in range(num_pathways)
        },
        index=adata.obs.index,
    )


def sim_gene_regulatory_networks(adata: sc.AnnData, outdir: str):
    """
    Generate two CSV files representing simulated GRNs:
    - sponge: uses Ensembl IDs (adata.var["ensemble_id"])
    - genie3: uses gene symbols (adata.var.index)
    """
    sponge_file = os.path.join(outdir, "sponge_network_filt.csv")
    genie_file = os.path.join(outdir, "genie_network_filt.csv")
    num_edges = 1000  # arbitrary number of edges per network

    # --- sponge network ---
    if "ensemble_id" in adata.var.columns:
        sponge_genes = adata.var["ensemble_id"].dropna().to_numpy()
    else:
        sponge_genes = np.array(
            [f"ENSG{100000+i}" for i in range(adata.shape[1])]
        )

    sponge_sources = np.random.choice(sponge_genes, size=num_edges)
    sponge_targets = np.random.choice(sponge_genes, size=num_edges)

    sponge_df = pd.DataFrame(
        {
            "geneA": sponge_sources,
            "geneB": sponge_targets,
            "df": 1,
            "cor": np.random.uniform(-0.3, 0.3, size=num_edges),
            "pcor": np.random.uniform(-0.3, 0.3, size=num_edges),
            "mscor": np.random.uniform(-0.1, 0.1, size=num_edges),
            "p.val": np.random.uniform(0.0, 1.0, size=num_edges),
            "p.adj": np.random.uniform(0.0, 1.0, size=num_edges),
        }
    )
    sponge_df.to_csv(sponge_file, index=False)

    # --- genie3 network ---
    genie_genes = adata.var.index.to_numpy()

    genie_sources = np.random.choice(genie_genes, size=num_edges)
    genie_targets = np.random.choice(genie_genes, size=num_edges)
    genie_weights = np.random.uniform(0.0, 1.0, size=num_edges)

    pd.DataFrame(
        {
            "regulatoryGene": genie_sources,
            "targetGene": genie_targets,
            "weight": genie_weights,
        }
    ).to_csv(genie_file, index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description=__doc__,
    )
    parser.add_argument("--adata", type=str, required=True)
    parser.add_argument("--output", type=str, required=True)
    args = parser.parse_args()

    adata = sc.read_h5ad(args.adata)

    funcs = [
        sim_genewise_scores,
        sim_leiden_scores,
        sim_Rscores,
        sim_liana_scores,
    ]
    for func in funcs:
        func(adata)

    sim_gene_regulatory_networks(adata, args.output)

    adata.write_h5ad(os.path.join(args.output, "adata.h5ad"))
