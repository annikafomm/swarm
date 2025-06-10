
import argparse
import os
import spatialdata as sd
import scanpy as sc
import squidpy as sq

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

parser.add_argument('-centrality_scores', action='store_true', help='Compute centrality scores')
parser.add_argument('-co_occurrence', action='store_true', help='Compute co-occurrence probability')
parser.add_argument('-nhood_enrichment', action='store_true', help='Compute neighborhood enrichment')

parser.add_argument('-cluster', type=str, default='leiden', help='Cluster key in adata.obs (default: leiden)')


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

# TODO: Make preprocessing user-configurable
if args.filter:
    print("Filtering ...")
    sc.pp.filter_cells(adata, min_counts=5)
    sc.pp.filter_genes(adata, min_cells=3)

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
        # TODO: Make the flavor configurable
        print("Computing Leiden clusters ...")
        sc.tl.leiden(adata, flavor="leidenalg", n_iterations=2, directed=False)  

# Compute centrality scores
if args.centrality_scores:
    print("Computing centrality scores ...")
    sq.gr.centrality_scores(adata, cluster_key=args.cluster)

# Compute co-occurrence probability
if args.co_occurrence:
    print("Computing co-occurrence probability ...")
    sq.gr.co_occurrence(adata, cluster_key=args.cluster)

# Compute neighborhood enrichment
if args.nhood_enrichment:
    print("Computing neighborhood enrichment ...")
    sq.gr.nhood_enrichment(adata, cluster_key=args.cluster)
   

# Compute Moran's I 
if args.moranI:
    print("Computing Moran's I ...")
    sq.gr.spatial_autocorr(adata, mode="moran")

# Compute Geary's C
if args.gearyC:
    print("Computing Geary's C ...")
    sq.gr.spatial_autocorr(adata, mode="geary")


# save AnnData object in file
print("Saving AnnData object ...")
if not args.output.endswith('.h5ad'):
    args.output += '.h5ad'
adata.write(args.output)
