import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning, module="xarray_schema")

import argparse
import os
import scanpy as sc
import squidpy as sq
import pandas as pd
import numpy as np
import json
import time
import anndata as ad

# Allow writing nullable strings in anndata
ad.settings.allow_write_nullable_strings = True

from calc_scores import log_message, format_runtime


def _ensure_spatial_neighbors(adata: ad.AnnData, logfile: str) -> bool:
    if "spatial" not in adata.obsm:
        log_message("Skipping regulatory autocorrelation: adata.obsm['spatial'] is missing.", logfile, 2)
        return False

    if "spatial_connectivities" in adata.obsp:
        return True

    t0 = time.time()
    sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)
    log_message(
        f"Spatial neighbors computed for regulatory autocorrelation in {format_runtime(t0)}",
        logfile,
        2,
    )
    return True


def _pick_stat_column(df: pd.DataFrame, preferred: str) -> str | None:
    if preferred in df.columns:
        return preferred

    numeric_cols = [c for c in df.columns if pd.api.types.is_numeric_dtype(df[c])]
    if numeric_cols:
        return numeric_cols[0]

    return None


def compute_regulatory_autocorr(adata: ad.AnnData, logfile: str) -> None:
    regulatory_keys = [
        key
        for key in adata.obsm.keys()
        if key.endswith("_genie3") or key.endswith("_sponge")
    ]

    if len(regulatory_keys) == 0:
        log_message("No regulatory score matrices found for autocorrelation.", logfile, 2)
        return

    if not _ensure_spatial_neighbors(adata, logfile):
        return

    spatial_connectivities = adata.obsp.get("spatial_connectivities", None)
    spatial_distances = adata.obsp.get("spatial_distances", None)

    if spatial_connectivities is None:
        log_message("Skipping regulatory autocorrelation: spatial connectivities missing.", logfile, 2)
        return

    log_message(
        f"Computing Moran's I / Geary's C for {len(regulatory_keys)} regulatory score matrices ...",
        logfile,
        2,
    )

    for key in regulatory_keys:
        try:
            score_df = adata.obsm[key]
            if not isinstance(score_df, pd.DataFrame):
                score_df = pd.DataFrame(score_df, index=adata.obs_names)

            if score_df.shape[1] == 0:
                log_message(f"Skipping {key}: matrix has no columns.", logfile, 4)
                continue

            tmp = ad.AnnData(X=score_df.to_numpy(dtype=np.float32))
            tmp.obs_names = adata.obs_names.copy()
            tmp.var_names = pd.Index(score_df.columns.astype(str))
            tmp.obsp["spatial_connectivities"] = spatial_connectivities.copy()
            if spatial_distances is not None:
                tmp.obsp["spatial_distances"] = spatial_distances.copy()

            t0 = time.time()
            sq.gr.spatial_autocorr(
                tmp,
                mode="moran",
                n_perms=None,
                corr_method="fdr_bh",
                show_progress_bar=False,
            )
            moran_df = tmp.uns.get("moranI", None)
            if isinstance(moran_df, pd.DataFrame):
                adata.uns[f"{key}_moranI"] = moran_df

            sq.gr.spatial_autocorr(
                tmp,
                mode="geary",
                n_perms=None,
                corr_method="fdr_bh",
                show_progress_bar=False,
            )
            geary_df = tmp.uns.get("gearyC", None)
            if isinstance(geary_df, pd.DataFrame):
                adata.uns[f"{key}_gearyC"] = geary_df

            log_message(f"Autocorrelation for {key} computed in {format_runtime(t0)}", logfile, 4)
        except Exception as e:
            log_message(f"Failed autocorrelation for {key}: {str(e)}", logfile, 4)


def combine_files(filename, description, args, logfile):
    t0 = time.time()
    adata_path = os.path.join(args.indir, filename)
    adata = sc.read_h5ad(adata_path)
    log_message(f"AnnData object loaded in {format_runtime(t0)}", logfile, 2)

    t0 = time.time()
    if description == "tg":
        scores_path = os.path.join(args.indir, "Rscores_tg")
    elif description == "st":
        scores_path = os.path.join(args.indir, "Rscores_st")

    if not os.path.isdir(scores_path):
        log_message(f"Scores directory not found: {scores_path}", logfile, 1)
        return

    merged_csv_keys = []
    merged_json_keys = []

    for filename in os.listdir(scores_path):
        if filename.endswith('.csv'):
            file_path = os.path.join(scores_path, filename)
            df = pd.read_csv(file_path, index_col=0)

            df_name = filename.replace(".csv", "")
            if df_name in adata.obsm.keys():
                log_message(f"The element with name {df_name} in obsm is overwritten.", logfile, 2)
            adata.obsm[df_name] = df.T
            merged_csv_keys.append(df_name)

        elif filename.endswith('.json'):
            file_path = os.path.join(scores_path, filename)
            with open(file_path, "r") as f:
                data_dict = json.load(f)

            df_name = filename.replace(".json", "")
            if df_name in adata.uns.keys():
                log_message(f"The element with name {df_name} in uns is overwritten.", logfile, 2)
            adata.uns[df_name] = data_dict
            merged_json_keys.append(df_name)

    if len(merged_csv_keys) == 0 and len(merged_json_keys) == 0:
        log_message(f"No R score files found in {scores_path}.", logfile, 1)
        return

    missing_csv = [k for k in merged_csv_keys if k not in adata.obsm.keys()]
    if len(missing_csv) > 0:
        raise RuntimeError(
            "R score merge failed: missing key(s) in adata.obsm after merge: "
            + ", ".join(missing_csv)
        )

    log_message(
        f"Merged {len(merged_csv_keys)} CSV and {len(merged_json_keys)} JSON score file(s) into AnnData.",
        logfile,
        2,
    )

    log_message(f"R score files loaded and added to the AnnData object in {format_runtime(t0)}", logfile, 2)

    print(adata)

    t0 = time.time()
    compute_regulatory_autocorr(adata, logfile)
    log_message(
        f"Regulatory autocorrelation post-processing finished in {format_runtime(t0)}",
        logfile,
        2,
    )

    # save AnnData object in file
    log_message("Saving AnnData object ...", logfile, 2)
    t0 = time.time()
    adata.write(adata_path)
    log_message(f"AnnData object written in {format_runtime(t0)}", logfile, 4)


def combine_files_multiome(filename, args, logfile, adata_map_path):
    scores_path = os.path.join(args.indir, "multiome")
    if not os.path.isdir(scores_path):
        log_message(f"Scores directory not found: {scores_path}", logfile, 1)
        return

    if not os.path.isfile(adata_map_path):
        log_message(f"Mapping file not found: {adata_map_path}. Skipping multiome score mapping.", logfile, 1)
        return

    t0 = time.time()
    adata_path = os.path.join(args.indir, filename)
    adata = sc.read_h5ad(adata_path)
    adata_map = sc.read_h5ad(adata_map_path)
    log_message(f"AnnData object loaded in {format_runtime(t0)}", logfile, 2)
    t0 = time.time()

    chromvar_found = False
    files_in_dir = os.listdir(scores_path)
    log_message(f"Files in multiome directory: {files_in_dir}", logfile, 2)

    for filename in files_in_dir:
        if filename.lower() == "peak_stats.csv":
            file_path = os.path.join(scores_path, filename)
            df = pd.read_csv(file_path)

            # clean column names and string values a bit
            df.columns = df.columns.str.strip()
            for col in df.select_dtypes(include="object").columns:
                df[col] = df[col].astype(str).str.strip()

            desired_cols = [
                "Gene",
                "Cluster",
                "Annotation",
                "Peak",
                "Class",
                "Link Score",
                "Link Z",
                "Link P",
                "Acc. T-stat",
                "Acc. FDR",
                "Accessible Cells",
                "P(expr|acc), cluster",
                "P(expr|acc), bg",
                "P(expr & acc), cluster",
                "P(expr & acc), all",
                "Enrichment, cluster",
                "Enrichment, all",
                "Delta P(expr|acc)",
                "Promoter Peaks",
                "Distal Peaks",
                "Pass Type",
            ]

            existing_cols = [col for col in desired_cols if col in df.columns]
            remaining_cols = [col for col in df.columns if col not in existing_cols]
            df = df[existing_cols + remaining_cols]

            if "peak_stats" in adata.uns:
                log_message(
                    "adata.uns['peak_stats'] is overwritten.",
                    logfile,
                    2
                )
            peak_stats_dict = df.to_dict(orient="list")
            adata.uns["peak_stats_dict"] = peak_stats_dict
            adata.uns["peak_stats"] = df
            log_message(
                "peak_stats added to adata.uns['peak_stats']",
                logfile,
                2
            )
        
        # motif stats
        if filename.lower() == "motif_stats.csv":
            file_path = os.path.join(scores_path, filename)
            df = pd.read_csv(file_path)

            # clean column names and string values
            df.columns = df.columns.str.strip()
            for col in df.select_dtypes(include="object").columns:
                df[col] = df[col].astype(str).str.strip()

            desired_cols = [
                "Gene",
                "Cluster",
                "TF",
                "Motif",
                "Prox Motif count",
                "Prox Bg count",
                "Prox Log2FC",
                "Prox p-value adj",
                "Dist Motif count",
                "Dist Bg count",
                "Dist Log2FC",
                "Dist p-value adj",
                "Prom Motif count",
                "FP Score",
                "Bg FP Score",
                "FP p-value adj",
                "Bg Size",
                "Flank sd",
                "Bg Flank sd",
                "Left Flank != 0",
                "Right Flank != 0",
            ]

            existing_cols = [col for col in desired_cols if col in df.columns]
            remaining_cols = [col for col in df.columns if col not in existing_cols]
            df = df[existing_cols + remaining_cols]

            for col in [
                "Cluster",
                "Prox Motif count",
                "Prox Bg count",
                "Prox Log2FC",
                "Prox p-value adj",
                "Dist Motif count",
                "Dist Bg count",
                "Dist Log2FC",
                "Dist p-value adj",
                "Prom Motif count",
                "FP Score",
                "Bg FP Score",
                "FP p-value adj",
                "Bg Size",
                "Flank sd",
                "Bg Flank sd",
                "Left Flank != 0",
                "Right Flank != 0",
            ]:
                if col in df.columns:
                    df[col] = pd.to_numeric(df[col], errors="coerce")

            if "motif_stats" in adata.uns:
                log_message(
                    "adata.uns['motif_stats'] is overwritten.",
                    logfile,
                    2
                )

            motif_stats_dict = df.to_dict(orient="list")
            adata.uns["motif_stats"] = df
            adata.uns["motif_stats_dict"] = motif_stats_dict
            log_message(
                "motif_stats added to adata.uns['motif_stats']",
                logfile,
                2
            )
        
        if filename.lower() == "chromvar_scores.csv":
            file_path = os.path.join(scores_path, filename)
            chromvar_scores = pd.read_csv(file_path, index_col=0)
            # adata is the AnnData object containing the mapping to spatial data (probabilities)
            M = adata_map.X # cells in dissociated scRNA assay x spatial cells probabilities
            S = np.array(chromvar_scores)[:,1:] # scores for dissociated data  cells x motifs

            M = np.asarray(M, dtype=np.float32)
            S = np.asarray(S, dtype=np.float32)

            spot_chromvar = M.T @ S # spatial cells x motifs
            motif_names = chromvar_scores.columns[1:]

            adata.obsm["chromvar_spot_scores"] = spot_chromvar
            adata.uns["chromvar_motifs"] = list(motif_names)
            chromvar_found = True
            log_message(f"chromvar_spot_scores added to adata.obsm", logfile, 2)

        # differential motif activity
        if filename.lower() == "differential_motif_activity.json":
            file_path = os.path.join(scores_path, filename)
            with open(file_path, "r") as f:
                data_dict = json.load(f)

            adata.uns["differential_motif_activity"] = data_dict
            log_message(f"differential_motif_activity added to adata.uns", logfile, 2)


        if filename.lower().startswith("diff_motif_activity_top_motifs_") and filename.lower().endswith(".csv"):
            file_path = os.path.join(scores_path, filename)
            df = pd.read_csv(file_path)

            desired_cols = ["motif", "p_val", "avg_diff", "pct.1", "pct.2", "p_val_adj"]
            existing_cols = [col for col in desired_cols if col in df.columns]
            remaining_cols = [col for col in df.columns if col not in existing_cols]
            df = df[existing_cols + remaining_cols]

            comparison = filename.replace(".csv", "").replace("diff_motif_activity_top_motifs_", "")

            if "diff_motif_activity_top_motifs" not in adata.uns:
                adata.uns["diff_motif_activity_top_motifs"] = {}

            if comparison in adata.uns["diff_motif_activity_top_motifs"]:
                log_message(
                    f"The comparison {comparison} in adata.uns['diff_motif_activity_top_motifs'] is overwritten.",
                    logfile,
                    2
                )

            adata.uns["diff_motif_activity_top_motifs"][comparison] = df
            log_message(
                f"{comparison} added to adata.uns['diff_motif_activity_top_motifs']",
                logfile,
                2
            )

    if not chromvar_found:
        log_message(f"Warning: chromvar_scores.csv not found in {scores_path}", logfile, 1)

    log_message(f"R score files loaded and added to the AnnData object in {format_runtime(t0)}", logfile, 2)

    # save AnnData object in file
    log_message("Saving AnnData object ...", logfile, 2)
    t0 = time.time()
    adata.write(adata_path)
    log_message(f"AnnData object written in {format_runtime(t0)}", logfile, 4)


def main():
    # Parser

    parser = argparse.ArgumentParser(description="Add R scores to adata file.")

    # input file paths
    parser.add_argument('-indir', type=str, required=True, help='Input dir file path')
    parser.add_argument('-log', type=str, required=True, help='Path to the log file')
    parser.add_argument('-Rscores', action='store_true', help='Whether to merge R scores or not')
    parser.add_argument('-multiome', action='store_true', help='Whether to merge multiome scores or not')

    args = parser.parse_args()

    logfile = args.log
    log_message(f"Combining scores pipeline started at {time.strftime('%Y-%m-%d %H:%M:%S')}", logfile)

    # Load the data
    log_message("Loading score data ...", logfile)

    adatas = sorted([
        filename
        for filename in os.listdir(args.indir)
        if filename.endswith('st_scores.h5ad') or filename.endswith('tg_scores.h5ad')
    ])

    if len(adatas) <= 0:
        log_message(f"The network scores could not be added to the AnnData object, because there are no .h5ad files in {args.indir}.", logfile)
    else:
        if args.multiome:
            adata_map_path = os.path.join(args.indir, "adata_map.h5ad")

        for filename in adatas:
            if filename.endswith('tg_scores.h5ad'):
                log_message("Combining files for the Tangram output ...", logfile)
                if args.Rscores:
                    combine_files(filename, "tg", args, logfile)
                if args.multiome:
                    combine_files_multiome(filename, args, logfile, adata_map_path)
            else:
                log_message("Combining files for the Spatial data ...", logfile)
                if args.Rscores:
                    combine_files(filename, "st", args, logfile)
                if args.multiome:
                    combine_files_multiome(filename, args, logfile, adata_map_path)

    log_message(f"Combining scores pipeline finished at {time.strftime('%Y-%m-%d %H:%M:%S')}\n", logfile)


if __name__ == "__main__":
    main()
