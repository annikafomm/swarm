import scanpy as sc
import tangram as tg

if __name__ == "__main__":
    ad_sp = sc.read_h5ad("/home/noah/Downloads/GSM6592052_M5.h5ad")
    ad_sc = sc.read("/home/noah/Downloads/Wu_annotated.h5ad")

    # ad_sp.var["gene_ids"] = ad_sp.var.index
    # ad_sp.var.set_index(ad_sp.var["feature_name"], drop=True, inplace=True)
    ad_sp = ad_sp[ad_sp.obs["in_tissue"] == 1]

    # Filter out mitochondrial genes
    ad_sp.var["MT_gene"] = [gene.startswith("MT-") for gene in ad_sp.var.index]
    ad_sp.obsm["MT"] = ad_sp[:, ad_sp.var["MT_gene"].values].X.toarray()
    ad_sp = ad_sp[:, ~ad_sp.var["MT_gene"].values]

    tg.pp_adatas(ad_sc, ad_sp, genes=None)

    ad_map = tg.map_cells_to_space(
        ad_sc, ad_sp, mode="clusters", cluster_label="cell_subclass"
    )

    ad_ge = tg.project_genes(ad_map, ad_sc, cluster_label="cell_subclass")

    tg.project_cell_annotations(ad_map, ad_ge, annotation="cell_subclass")
    ad_ge.obsm["spatial"] = ad_sp.obsm["spatial"]
    ad_ge.uns["spatial"] = ad_sp.uns["spatial"]

    ad_ge.write("/home/noah/Downloads/tangram_adata.h5ad")
