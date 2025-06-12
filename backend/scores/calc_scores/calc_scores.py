import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

import argparse
import os
import spatialdata as sd
import scanpy as sc
import squidpy as sq


def main():
    # Parser

    parser = argparse.ArgumentParser(description="Calculate spatial scores for spatial omics data.")

    # input and output file paths
    parser.add_argument('-input', type=str, required=True, help='Input AnnData file path')
    parser.add_argument('-output', type=str, default='output.h5ad', help='Output AnnData file path')

    # preprocessing options
    parser.add_argument('-filter', action='store_true', help='Apply filtering')
    parser.add_argument('-normalize', action='store_true', help='Apply normalization')

    # possible scores
    parser.add_argument('-moranI', action='store_true', help='Compute Moran\'s I')
    parser.add_argument('-gearyC', action='store_true', help='Compute Geary\'s C')

    parser.add_argument('-centrality_scores', action='store_true', help='Compute all centrality scores (closeness, average clustering, degree)')    
    parser.add_argument('-co_occurrence', action='store_true', help='Compute co-occurrence probability')
    parser.add_argument('-nhood_enrichment', action='store_true', help='Compute neighborhood enrichment')

    # additional options
    parser.add_argument('-cluster', type=str, default='leiden', help='Cluster key in adata.obs (default: leiden)')
    parser.add_argument('-n_perms', type=int, default=None, help='Number of permutations for scores')
    parser.add_argument('-n_jobs', type=int, default=None, help='Number of jobs for parallel processing')

    #parser.add_argument('-corr_method', type=str, default='fdr_bh', help='Correlation method for Moran\'s I and Geary\'s C (default: benjamini-hochberg)')
    #parser.add_argument('-numba_parallel', action='store_true', help='Use numba for parallel processing')
    #parser.add_argument('-interval', type=int, default=50, help='Distances interval at which co-occurrence is computed (default: 50)')
    #parser.add_argument('-n_splits', type=int, default=None, help='Number of splits in which to divide the spatial coordinates')

    args = parser.parse_args()


    # Load the data
    print("Reading input file ...")
    if not os.path.exists(args.input):
        raise FileNotFoundError(f"Input file {args.input} does not exist.")

    if args.input.endswith('.h5ad'):
        adata = sc.read_h5ad(args.input)
    elif args.input.endswith('.zarr'):
        sdata = sd.read_zarr(args.input)
        adata = sdata.tables["table"]
    else:
        raise ValueError("Unsupported file format. Please provide a .h5ad or .zarr file.")


    # Preprocessing

    if args.filter:
        print("Filtering ...")
        sc.pp.filter_cells(adata, min_counts=10) # GitHub CoPilot: For spatial transcriptomics (like Xenium) typical values are 10–100.
        # GitHub CoPilot: For single cell data typical values are 200–500.
        sc.pp.filter_genes(adata, min_cells=3) # GitHub CoPilot: For spatial data (like Xenium), 3–10 is typical.

    if args.normalize:
        print("Normalization ...")
        adata.layers["counts"] = adata.X.copy()

        sc.pp.normalize_total(adata, inplace=True) # Normalize counts per cell
        sc.pp.log1p(adata) # Logarithmize
        sc.pp.pca(adata) # do principal component analysis


    # Calculate spatial scores

    if 'spatial_connectivities' not in adata.obsp.keys() or 'spatial_distances' not in adata.obsp.keys():
        print("Computing the spatial neighbors ...")
        sq.gr.spatial_neighbors(adata, coord_type="generic", delaunay=True)


    # check if the cluster key exists in adata.obs if needed
    if args.centrality_scores or args.co_occurrence or args.nhood_enrichment:
        if args.cluster not in adata.obs.keys() and args.cluster != "leiden":
            raise ValueError(f"Cluster key '{args.cluster}' not found in adata.obs. Please provide a valid cluster key.")
        
        # Compute Leiden clusters if not already present
        if args.cluster == "leiden" and "leiden" not in adata.obs.keys():
            print("Computing the neighborhood graph (leiden cluster) ...")
            sc.pp.neighbors(adata) # compute a neighborhood graph of the observations
            print("Computing Leiden clusters ...")
            sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False) # not user configurable, because makeshift solution for when no cluster key is provided

    # Compute centrality scores
    if args.centrality_scores:
        print("Computing centrality scores ...")
        sq.gr.centrality_scores(adata, cluster_key=args.cluster, score=None, n_jobs=args.n_jobs, show_progress_bar=True)

    # Compute co-occurrence probability
    if args.co_occurrence:
        print("Computing co-occurrence probability ...")
        sq.gr.co_occurrence(adata, cluster_key=args.cluster, n_jobs=args.n_jobs, show_progress_bar=True)

    # Compute neighborhood enrichment
    if args.nhood_enrichment:
        print("Computing neighborhood enrichment ...")
        sq.gr.nhood_enrichment(adata, cluster_key=args.cluster, seed=0, n_perms=3 if args.n_perms is None or args.n_perms < 3 else args.n_perms, n_jobs=args.n_jobs, show_progress_bar=True)
    

    # Compute Moran's I 
    if args.moranI:
        print("Computing Moran's I ...")
        sq.gr.spatial_autocorr(adata, mode="moran", seed=0, n_perms=args.n_perms, n_jobs=args.n_jobs, transformation=args.n_perms is None, show_progress_bar=True)

    # Compute Geary's C
    if args.gearyC:
        print("Computing Geary's C ...")
        sq.gr.spatial_autocorr(adata, mode="geary", seed=0, n_perms=args.n_perms, n_jobs=args.n_jobs, transformation=args.n_perms is None, show_progress_bar=True)


    # save AnnData object in file
    print("Saving AnnData object ...")
    if not args.output.endswith('.h5ad'):
        args.output += '.h5ad'
    adata.write(args.output)

if __name__ == "__main__":
    main()
    print("Done")