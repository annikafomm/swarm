import scanpy as sc
import time

#adata = sc.read_h5ad(r"C:\RebEll\Computer\BioinformatikMaster\Semester8\MaPra\calc_scores\xenium.h5ad")
adata = sc.read_h5ad(r"C:\RebEll\Computer\BioinformatikMaster\Semester8\MaPra\rebecca\datasets\GSM6592055_M8.h5ad")

sc.pp.filter_cells(adata, min_counts=10) 
sc.pp.filter_genes(adata, min_cells=5)

sc.pp.normalize_total(adata, inplace=True) # Normalize counts per cell
sc.pp.log1p(adata) # Logarithmize
sc.pp.pca(adata) # do principal component analysis

sc.pp.neighbors(adata)


start = time.time()
"""
scanpy.tl.leiden(adata, ...)
- resolution=1             --> higher values lead to more clusters
- random_state=0           --> use an integer for reproducibility
- n_iterations=2           --> -1 runs until convergence; set to a positive integer to limit iterations
- directed=False           --> True for directed graphs, False for undirected graphs
- flavor='igraph'          --> 'igraph' or 'leidenalg'
"""

sc.tl.leiden(adata, flavor="igraph", n_iterations=2, directed=False) # took 32.13 seconds
#sc.tl.leiden(adata, flavor="igraph", resolution=5, n_iterations=2, directed=False) # took 29.97 seconds
#sc.tl.leiden(adata, flavor="igraph", resolution=20, n_iterations=2, directed=False) # took 31.67 seconds

#sc.tl.leiden(adata, flavor="leidenalg", n_iterations=2, directed=False) # took 0.24 seconds
#sc.tl.leiden(adata, flavor="leidenalg", resolution=5, n_iterations=2, directed=False) # took 0.22 seconds
#sc.tl.leiden(adata, flavor="leidenalg", resolution=20, n_iterations=2, directed=False) # took 0.26 seconds

end = time.time()
print(f"Leiden clustering took {end - start:.2f} seconds")
