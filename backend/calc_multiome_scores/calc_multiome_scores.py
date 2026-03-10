import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="xarray_schema")

import argparse
import os
import scanpy as sc
import squidpy as sq
from scipy import io, sparse
import pandas as pd
import time
import anndata as ad

# Allow writing nullable strings in anndata
ad.settings.allow_write_nullable_strings = True

import sys
sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

# python3 calc_multiome_scores.py --dir /workspaces/swarm/backend/uploads/job_1765267595498_f19ae350-7f14-4fdf-88d7-77a40d81c2c4/plasmidpoop --log /workspaces/swarm/backend/uploads/multiome.txt -chromvar -moranI


def log_message(msg, logfile, indent=0):
    prefix = " " * indent
    line = f"{prefix}{msg}"
    print(line)  # console
    with open(logfile, "a") as f:
        f.write(line + "\n")  # write once, no extra blank lines

def format_runtime(t0):
    elapsed = time.time() - t0  # seconds since start
    mins = int(elapsed // 60)
    secs = int(round(elapsed % 60))
    return f"{mins} min {secs} sec"


def _ensure_spatial_neighbors(motif_adata, base_adata=None,
                              coord_type="generic"):
    """
    Make sure motif_adata has spatial coords + neighbor graph.
    - If base_adata is given, copy what we can from there.
    - If any pieces are missing, compute them with squidpy.
    """
    # copy coordinates from base_adata if possible
    if base_adata is not None and "spatial" in base_adata.obsm:
        motif_adata.obsm["spatial"] = base_adata.obsm["spatial"].copy()

    # copy neighbor graph pieces if available
    if base_adata is not None:
        for key in ("spatial_connectivities", "spatial_distances"):
            if key in base_adata.obsp and key not in motif_adata.obsp:
                motif_adata.obsp[key] = base_adata.obsp[key].copy()

        if "spatial_neighbors" in base_adata.uns and "spatial_neighbors" not in motif_adata.uns:
            motif_adata.uns["spatial_neighbors"] = base_adata.uns["spatial_neighbors"].copy()

    # if anything still missing, compute with squidpy
    need_connectivities = any(
        key not in motif_adata.obsp
        for key in ("spatial_connectivities", "spatial_distances")
    )
    need_uns = "spatial_neighbors" not in motif_adata.uns

    if need_connectivities or need_uns:
        print("Computing spatial neighbors...")
        if "spatial" not in motif_adata.obsm:
            raise ValueError("No .obsm['spatial'] found to compute spatial neighbors from.")
        sq.gr.spatial_neighbors(
            motif_adata,
            coord_type=coord_type   # "grid" for Visium, "generic" for others
        )


def spatial_correlation(adata, args, logfile):

    # create new adata for this:
    # motif activity: (n_cells, n_motifs)
    motif_mat = adata.obsm["chromvar_spot_scores"].copy()
    motif_names = adata.uns["chromvar_motifs"].copy()  # should be length n_motifs

    motif_adata = ad.AnnData(
        X=motif_mat,
        obs=adata.obs.copy(),
        var=pd.DataFrame(index=motif_names),
    )

    _ensure_spatial_neighbors(
        motif_adata,
        base_adata=adata,
        coord_type="generic"   # or "generic" depending on your data
    )

    # compute spatial correlation
    # TODO: DOES IT MAKE SENSE THAT THE MOTIFS OCCUR MULTIPLE TIMES IF THEY MAP TO MULTIPLE TFs?
    if args.moranI or args.gearyC:
        motif_map = pd.read_csv(os.path.join(args.dir, "motif_to_tf.csv"))
        motif_map_merit_exp = (
            motif_map
            .assign(TF=motif_map["TF"].str.split("::"))
            .explode("TF")
            .dropna(subset=["TF"])
        )

    # Compute Moran's I
    if args.moranI:
        print("Computing Moran's I...")
        t0 = time.time()
        sq.gr.spatial_autocorr(motif_adata,
                                mode="moran",
                                seed=42,
                                n_perms=args.n_perms_autocorr_mI,
                                transformation=args.n_perms_autocorr_mI is None,
                                two_tailed = args.two_tailed_mI,
                                corr_method = args.corr_method_mI,
                                show_progress_bar=True)
        log_message(f"Moran's I scores calculated in {format_runtime(t0)}", logfile, 2)

        # map tf to motif
        merged_table = motif_adata.uns["moranI"].merge(motif_map_merit_exp, left_index=True, right_on="motif_id")
        # order columns
        merged_table = merged_table[["TF", "motif_id", "I", "pval_norm", "var_norm", "pval_z_sim", "pval_sim", "var_sim", "pval_norm_fdr_bh","pval_z_sim_fdr_bh", "pval_sim_fdr_bh"]]
        merged_table = merged_table.sort_values(by="TF")
        # store in adata.uns
        adata.uns["chromvar_moranI"] = merged_table


    # Compute Geary's C
    if args.gearyC:
        t0 = time.time()
        sq.gr.spatial_autocorr(motif_adata,
                                mode="geary",
                                seed=42,
                                n_perms=args.n_perms_autocorr_gC,
                                transformation=args.n_perms_autocorr_gC is None,
                                two_tailed = args.two_tailed_gC,
                                corr_method = args.corr_method_gC,
                                show_progress_bar=True)
        log_message(f"Geary's C scores calculated in {format_runtime(t0)}", logfile, 2)

        # map tf to motif
        merged_table = motif_adata.uns["gearyC"].merge(motif_map_merit_exp, left_index=True, right_on="motif_id")
        # order columns
        merged_table = merged_table[["TF", "motif_id", "C", "pval_norm", "var_norm", "pval_z_sim", "pval_sim", "var_sim", "pval_norm_fdr_bh","pval_z_sim_fdr_bh", "pval_sim_fdr_bh"]]
        merged_table = merged_table.sort_values(by="TF")
        # store in adata.uns
        adata.uns["chromvar_gearyC"] = merged_table


    t0 = time.time()
    #filename = "adata_tg_scores_w_chromvar.h5ad"
    filename = "adata_tg_scores.h5ad"
    adata.write(os.path.join(args.dir, filename))
    filename = "adata_st_scores.h5ad"
    adata.write(os.path.join(args.dir, filename))
    log_message(f"AnnData object written in {format_runtime(t0)}", logfile, 4)


def main():
    # Parser
    parser = argparse.ArgumentParser(description="Calculate spatial scores for spatial omics data.")

    # input and output file paths
    parser.add_argument('--dir', type=str, required=True, help='Input AnnData file path')
    parser.add_argument('--log', type=str, required=True, help='Path to the log file')

    parser.add_argument('-chromvar', action='store_true', help='contains chromVar scores')

    parser.add_argument('-moranI', action='store_true', help='Compute Moran\'s I')
    parser.add_argument('-two_tailed_mI', action='store_true', help='two-tailed test for Moran\'s I')
    parser.add_argument('-corr_method_mI', type=str, default="fdr_bh", help='correction method for Moran\'s I')
    parser.add_argument('-n_perms_autocorr_mI', type=int, default=1000, help='')

    parser.add_argument('-gearyC', action='store_true', help='Compute Geary\'s C')
    parser.add_argument('-two_tailed_gC', action='store_true', help='two-tailed test for Geary\'s C')
    parser.add_argument('-corr_method_gC', type=str, default="fdr_bh", help='correction method for Geary\'s C')
    parser.add_argument('-n_perms_autocorr_gC', type=int, default=1000, help='')

    parser.add_argument('-differential_motif_activity', action='store_true', help='contains differential motif activity scores')
    parser.add_argument('-motif_enrichment', action='store_true', help='contains motif enrichment scores')
    parser.add_argument('-footprinting', action='store_true', help='contains footprinting scores')

    args = parser.parse_args()

    adata = sc.read_h5ad(os.path.join(args.dir, "adata_tg_scores.h5ad"))

    if args.chromvar:
        spatial_correlation(adata, args, args.log)



if __name__ == "__main__":
    main()
