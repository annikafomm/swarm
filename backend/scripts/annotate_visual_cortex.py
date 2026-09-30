"""
One-off enrichment for the Visual Cortex (Visium) builtin dataset:

1. Cell type: majority-vote the existing Stereoscope-derived `obsm['proportions_class']`
   composition estimates into a single `obs['cell_type']` categorical call per spot, and also
   write each class's raw proportion as its own `obs[<class>]` column (matching the convention
   already used by the Heart dataset's deconvolution fractions), so both the vote and the
   underlying evidence are visible.
2. Region: since no real anatomical layer/region annotation exists for this dataset (the only
   one found upstream, allen-1-meta.tsv, is a hand-drawn binary isocortex mask the original
   author himself calls "rather sloppy"), infer cortical layer per spot via scanpy's
   `score_genes` module-scoring against canonical mouse cortex layer marker panels, then argmax.
   This is a marker-based inference, not ground truth -- written to `obs['region']` with the
   per-layer scores kept as `obs['region_score_<layer>']` for transparency/QC.

Writes an enriched copy of the h5ad; does not mutate the original file in place.
"""

import sys
import scanpy as sc
import pandas as pd
import numpy as np

SRC = "/workspaces/swarm/backend/data/visual_cortex/dumpase1/adata_st_scores.h5ad"
DST = "/workspaces/swarm/backend/data/visual_cortex/dumpase1/adata_st_scores_annotated.h5ad"

# Canonical mouse cortex layer marker panels (well-established literature markers; Allen Institute
# ISH atlas / Tasic et al. 2018 nomenclature). Non-cortical/white-matter markers included so spots
# outside the cortical column don't get forced into a cortical layer call.
LAYER_MARKERS = {
    "L1": ["Reln", "Ndnf", "Lamp5"],
    "L2/3": ["Cux2", "Calb1", "Rasgrf2"],
    "L4": ["Rorb", "Rspo1", "Scnn1a"],
    "L5": ["Fezf2", "Deptor", "Bcl11b", "Etv1"],
    "L6": ["Foxp2", "Tle4", "Ntsr1", "Ctgf"],
    "White matter": ["Mbp", "Mobp", "Plp1"],
}


def main():
    print(f"Loading {SRC} ...", flush=True)
    adata = sc.read_h5ad(SRC)
    print("Loaded:", adata.shape, flush=True)

    if "proportions_class" not in adata.obsm:
        sys.exit("obsm['proportions_class'] not found -- aborting.")

    props = adata.obsm["proportions_class"]
    if not isinstance(props, pd.DataFrame):
        sys.exit("obsm['proportions_class'] is not a DataFrame -- can't recover class names, aborting.")

    print("proportions_class columns:", list(props.columns), flush=True)

    # ----- 1. Cell type: majority vote + raw proportions -----
    collisions = [c for c in props.columns if c in adata.obs.columns]
    if collisions:
        sys.exit(f"obs already has columns colliding with class names: {collisions} -- aborting, resolve manually.")

    adata.obs["cell_type"] = props.idxmax(axis=1).astype(str).values
    for col in props.columns:
        adata.obs[col] = props[col].values
    print("Added obs['cell_type'] (majority vote) and", len(props.columns), "proportion columns.", flush=True)
    print("cell_type value counts:\n", adata.obs["cell_type"].value_counts(), flush=True)

    # ----- 2. Region: marker-based cortical layer inference -----
    # Match markers case-insensitively and index with the dataset's own casing. This file's
    # symbols were upper-cased ("CUX2") until scripts/restore_var_symbols.py put the mm10 casing
    # back ("Cux2"), and LAYER_MARKERS below is written in mouse casing; going through
    # var_by_upper keeps this script working against either form.
    var_by_upper = {v.upper(): v for v in adata.var_names}
    score_cols = []
    for layer, genes in LAYER_MARKERS.items():
        present = [var_by_upper[g.upper()] for g in genes if g.upper() in var_by_upper]
        missing = [g for g in genes if g.upper() not in var_by_upper]
        if missing:
            print(f"[{layer}] missing from var_names, skipping: {missing}", flush=True)
        if not present:
            print(f"[{layer}] no markers present at all -- skipping this layer entirely.", flush=True)
            continue
        score_name = f"region_score_{layer.replace('/', '')}"
        sc.tl.score_genes(adata, present, score_name=score_name)
        score_cols.append((layer, score_name))
        print(f"[{layer}] scored using {present}", flush=True)

    if not score_cols:
        print("WARNING: no layer marker genes found -- skipping region inference, keeping cell_type only.", flush=True)
    else:
        score_df = pd.DataFrame({layer: adata.obs[col] for layer, col in score_cols})
        # score_genes' raw output isn't comparable across different marker panels -- a panel of
        # generally highly-expressed genes (e.g. myelin markers Mbp/Mobp/Plp1) can win the argmax
        # everywhere purely on scale, not true relative enrichment. Z-score each layer's score
        # across spots first so the comparison reflects each panel's own relative enrichment.
        zscored = (score_df - score_df.mean()) / score_df.std(ddof=0)
        adata.obs["region"] = zscored.idxmax(axis=1).values
        print("region value counts:\n", adata.obs["region"].value_counts(), flush=True)

    print(f"Writing {DST} ...", flush=True)
    adata.write_h5ad(DST)
    print("Done.", flush=True)


if __name__ == "__main__":
    main()
