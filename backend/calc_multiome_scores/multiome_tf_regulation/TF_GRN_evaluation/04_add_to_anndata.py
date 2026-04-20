
import os

import scanpy as sc
import pandas as pd
import anndata as ad

adata_tg = "path/to/adata_tg.h5ad"
outpath = "path/to/output_directory"

adata = sc.read_h5ad(adata_tg)

file_path = os.path.join(outpath, "peak_stats.csv")
df = pd.read_csv(file_path)

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
peak_stats_dict = df.to_dict(orient="list")
key = "collagen_fibril_orga"
if key not in adata.uns or not isinstance(adata.uns[key], dict):
    adata.uns[key] = {}

adata.uns[key]["peak_stats"] = df
adata.uns[key]["peak_stats_dict"] = peak_stats_dict

file_path = os.path.join(outpath, "motif_stats.csv")
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

motif_stats_dict = df.to_dict(orient="list")

adata.uns["collagen_fibril_orga"]["motif_stats"] = df
adata.uns["collagen_fibril_orga"]["motif_stats_dict"] = motif_stats_dict


adata.write_h5ad(adata_tg)

    