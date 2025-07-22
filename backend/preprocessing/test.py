import scanpy as sc
import matplotlib
matplotlib.use("Agg")

adata = sc.datasets.pbmc3k()
sc.pp.calculate_qc_metrics(adata, inplace=True)
sc.pp.highly_variable_genes(adata, flavor='seurat_v3', n_top_genes=2000)
sc.pl.highly_variable_genes(adata, save="_test.png")
