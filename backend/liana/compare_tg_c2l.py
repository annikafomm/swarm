import matplotlib.pyplot as plt
import scanpy as sc
import squidpy as sq

if __name__ == "__main__":
    tg_ad = sc.read_h5ad("/home/noah/Downloads/liana_tangram_adata.h5ad")
    c2l_ad = sc.read_h5ad(
        "/home/noah/Downloads/liana_cell2location_adata.h5ad"
    )

    comp_tf_rels = [
        "Myoepithelial<->FOXF2",
        "NKT cells<->IRX5",
        "Plasmablasts<->ZNF440",
        # "DCs<->GLI2",
    ]

    sq.pl.spatial_scatter(
        tg_ad, color=comp_tf_rels, library_id="GSM6592052_M5", cmap="coolwarm"
    )
    plt.show()

    sq.pl.spatial_scatter(
        c2l_ad, color=comp_tf_rels, library_id="GSM6592052_M5", cmap="coolwarm"
    )
    plt.show()

    lr_cos_sim = "TNFSF9^TNFRSF9_cosine-similarity"
    lr_pval = "TNFSF9^TNFRSF9_p-value"
    lr_cat = "TNFSF9^TNFRSF9_category"

    sq.pl.spatial_scatter(
        c2l_ad, color=lr_cos_sim, library_id="GSM6592052_M5", cmap="magma"
    )
    plt.show()

    sq.pl.spatial_scatter(
        c2l_ad, color=lr_pval, library_id="GSM6592052_M5", cmap="magma_r"
    )
    plt.show()

    sq.pl.spatial_scatter(
        c2l_ad, color=lr_cat, library_id="GSM6592052_M5", cmap="coolwarm"
    )
    plt.show()

    nmf_factors = [f"NMF_Factor_{i}" for i in range(1, 4)]

    sq.pl.spatial_scatter(
        c2l_ad, color=nmf_factors, library_id="GSM6592052_M5"
    )
    plt.show()
