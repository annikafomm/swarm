import argparse
import os
import scanpy as sc
from gridding_mapping import expand_spot_adata_to_cells


def main():
    # Parser
    parser = argparse.ArgumentParser(description="Calculate back mapping for Xenium.")

    # input and output file paths
    parser.add_argument('-indir', type=str, required=True, help='Input dir file path')
    parser.add_argument('-log', type=str, required=True, help='Path to the log file')
    parser.add_argument('-tangram', action="store_true", help='Whether to use Tangram anndata for backmapping (default: False)')

    args = parser.parse_args()
    # Xenium only: scores were computed on the grid-level (adata_work),
    # so we broadcast grid-level obs/obsm back to each cell using the stored mapping.
    print(args.tangram)
    if args.tangram:
        work_path = os.path.join(args.indir, "adata_tg_scores.h5ad")
        adata_work = sc.read_h5ad(work_path)
    else:
        work_path = os.path.join(args.indir, "adata_st_scores.h5ad")
        adata_work = sc.read_h5ad(work_path)
    adata_cells_path = os.path.join(args.indir, "xenium_map.h5ad")
    adata_cells = sc.read_h5ad(adata_cells_path)
    print(adata_work)
    print(adata_cells)

    adata_cells = expand_spot_adata_to_cells(
        adata_scores=adata_work,
        adata_map=adata_cells,
        spot_col="grid_spot"
    )
    print(f"Broadcasted grid-level data to cells: {adata_cells} cells")

    adata_cells.write(
        os.path.join(args.indir, "xenium_cells_with_grid_scores.h5ad")
    )

if __name__ == "__main__":
    main()
