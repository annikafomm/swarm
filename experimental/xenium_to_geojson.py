# coding: utf-8
import geopandas
import matplotlib.pyplot as plt
import nichepca as npc
import scanpy as sc
import spatialdata as sd
import squidpy as sq
from scipy.stats import mode
from spatialdata_io import xenium

if __name__ == "__main__":
    # xenium_path = "/nfs/data3/mopitas/mapra/datasets/output-XETG00125__0042682__3802__20240426__200549"
    xenium_path = "/nfs/data3/mopitas/mapra/datasets/Janesick_outs/"
    sdata = xenium(xenium_path)
    adata_raw = sdata["table"]
    adata = adata_raw.copy()
    adata.obsm["spatial"]

    adata.obs["z_level"].value_counts()
    sq.pl.spatial_scatter(
        adata,
        shape="hex",
        color=["z_level"],
        img=False,
        save="spatial_scatter_zlvl_hex.png",
    )

    # Default parameters from tutorial
    npc.wf.nichepca(adata, knn=25)
    sc.pp.neighbors(adata, use_rep="X_npca")
    sc.tl.leiden(adata)
    adata

    adata.obs["leiden"].value_counts()
    sq.pl.spatial_scatter(
        adata,
        shape=None,
        color=["leiden"],
        legend_loc=None,
        save="spatial_scatter_npca_knn25_res01.png",
    )

    # Alternative parameters results in usable clustering
    # Higher resolution => more clusters
    adata = adata_raw.copy()
    knn = 25
    res = 100  # 100
    npc.wf.nichepca(adata, knn=knn)
    sc.pp.neighbors(adata, use_rep="X_npca")
    # sc.tl.leiden(adata, resolution=res, adjacency=adata.obsp["spatial_connectivities"])

    # Smoothing based on spatial graph, idea taken from
    # https://www.sc-best-practices.org/spatial/domains.html
    nn_graph_genes = adata.obsp["connectivities"]
    sq.gr.spatial_neighbors(adata)
    nn_graph_space = adata.obsp["spatial_connectivities"]
    alpha = 0.8
    joint_graph = (1 - alpha) * nn_graph_genes + alpha * nn_graph_space
    sc.tl.leiden(adata, resolution=res, adjacency=joint_graph)

    adata.obs["leiden"].value_counts()
    sq.pl.spatial_scatter(
        adata,
        shape=None,
        color=["leiden"],
        legend_loc=None,
        save=f"spatial_scatter_npca_knn{knn}_res{res}_alpha{alpha}.png",
    )

    # Test the nichePCA-based binning
    cell_boundaries = sdata.shapes["cell_boundaries"]
    gdf = adata.obs
    gdf.index = gdf["cell_id"]
    gdf = geopandas.GeoDataFrame(data=gdf, geometry=cell_boundaries["geometry"])

    # HACK: Buffer usually extends the size of the polygon by `distance`.
    # We use it here with distance 0 to avoid invalid geometries.
    # https://stackoverflow.com/questions/74779301/unable-to-assign-free-hole-to-a-shell-error-when-flattening-polygons
    gdf["geometry"] = gdf["geometry"].buffer(0)

    # Dissolve works like group_by but unions the polygons.
    # `method=coverage` for non-overlapping polygons, also is faster.
    gdf = gdf.dissolve(by="leiden", method="coverage")

    # The convex hull creates large polygons containing all inner polygons.
    gdf["geometry"] = gdf["geometry"].convex_hull

    gdf["tmp"] = gdf.index
    gdf.plot(column="tmp")
    plt.show()
    plt.savefig("/nfs/home/students/n.nussbaumer/test_leiden_smooth.png")
    plt.savefig("/nfs/home/students/n.nussbaumer/test_leiden_smooth.svg")

    gdf.to_file(f"xenium_alpha{alpha}.geojson")

    # gdf = gdf.dissolve(by="leiden", method="coverage")
    # gdf["geometry"] = gdf["geometry"].buffer(0.1).simplify(0.2)
    # gdf = gdf.union_all(method="unary").polygonize()
    # gdf["tmp"] = gdf.index
    # gdf.plot(column="tmp")
    # plt.show()
    #
    # # get list of geometries
    # geoms = gdf["geometry"].tolist()
    # # iterate over all combinations of polygons and get the intersections (overlaps)
    # overlaps = geopandas.GeoDataFrame(
    #     geopandas.GeoSeries(
    #         [
    #             poly[0].intersection(poly[1])
    #             for poly in itertools.combinations(geoms, 2)
    #             if poly[0].intersects(poly[1])
    #         ]
    #     ),
    #     columns=["geometry"],
    # )
    # # erase the overlaps from the original geodataframe
    # overlaps = overlaps.polygonize().to_frame()
    # gdf2 = gdf.overlay(overlaps, how="difference")
    #
    # gdf2.to_file("xenium_spatialn_nool.geojson")
    # gdf2.plot(column="tmp")
    # plt.show()
