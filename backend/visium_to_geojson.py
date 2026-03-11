import argparse
import os
import geopandas as gpd
import numpy as np
import pandas as pd
import scanpy as sc
import simplejson as json
from scipy import sparse

# python ../backend/visium_to_geojson.py --adata ../backend/uploads/job_1764657831787_merit/junkDNA420/adata_st_scores.h5ad --outpath /workspaces/swarm/backend/uploads/job_1764657831787_merit/junkDNA420/hexagons.geojson

grn_score_names = (
    "aucell_scores",
    "spongeffects_GSVA_scores",
    "spongeffects_ssGSEA_scores",
    "viper_scores",
    "grid_viper_scores",
    "grid_aucell_scores",
    "grid_spongeffects_GSVA_scores",
    "grid_spongeffects_ssGSEA_scores",
)

genie3_score_names = [name + "_genie3" for name in grn_score_names]
sponge_score_names = [
    name + "_sponge"
    for name in grn_score_names
    if not name.startswith("viper")  # Viper score not available for sponge
]

genewise_scores = ["moranI", "gearyC"]
motifwise_scores = ["chromvar_moranI", "chromvar_gearyC"]
cluster_wise_scores = ["diff_motif_activity_top_motifs"]

class Hexagons:
    """
    Class to convert Visium or Transformed Xenium spatial data to GeoJSON format with hexagonal geometry.
    This class handles the parsing of AnnData objects, Spatial scores, Gene regulatory Scores and spatial coordinates
    to create a GeoJSON representation of the spatial data with annotated properties.
    """

    def __init__(
        self,
        anndata,
        radius=5,
        scale=0.1,
        data_type="visium",
        motif_groups=None,
    ):

        self.anndata = anndata
        self.radius = radius
        self.type = "FeatureCollection"
        self.geometry_type = "Polygon"
        self.scale = scale
        self.data_type = data_type

        self.global_scores = {}
        self.motif_groups = motif_groups or {}
        self.motif_names = list(anndata.uns["chromvar_motifs"]) if "chromvar_motifs" in anndata.uns else []
        # mapping for case-insensitive gene lookups
        self.var_upper_to_name = {
            g.upper(): g for g in self.anndata.var_names
        }
        self.coordinates, self.centers = self.parse_coordinates()
        self.obs = self.parse_obs()

    def parse_obs(self):
        obs = self.anndata.obs
        properties = obs.copy()
        # Convert to dict with id as key and properties as values
        properties_dict = properties.to_dict(orient="index")
        return properties_dict

    def hexagon_points(self, x, y, radius):
        # Scale centers
        x = x * self.scale
        y = y * self.scale
        return [
            (
                x + radius * np.cos(np.pi / 3 * i),
                y + radius * np.sin(np.pi / 3 * i),
            )
            for i in range(6)
        ] + [(x + radius * np.cos(0), y + radius * np.sin(0))]

    def parse_coordinates(self):
        anndata_spatial_coordinates = self.anndata.obsm["spatial"].copy()
        barcodes = self.anndata.obs.index

        hex_coords = []
        coords = []
        if self.data_type == "visium":
            for coord_tuple in anndata_spatial_coordinates:
                x, y = coord_tuple
                hexagon = self.hexagon_points(x=x, y=y, radius=self.radius)
                hex_coords.append(hexagon)
                coords.append([int(x) * self.scale, int(y) * self.scale])
            barcodes = self.anndata.obs.index

        # Xenium mode:
        # Render each cell as a small hexagon instead of a large spot.
        # Grid-level scores are already broadcasted into obs/obsm and will be attached
        # to each cell feature for visualization.

        elif self.data_type == "xenium":
            small_r = getattr(self, "cell_radius", None)
            if small_r is None:
                small_r = self.radius * 0.2

            for coord_tuple in anndata_spatial_coordinates:
                x, y = coord_tuple
                hexagon = self.hexagon_points(x=x, y=y, radius=small_r)
                hex_coords.append(hexagon)
                coords.append([float(x) * self.scale, float(y) * self.scale])
        else:
            raise ValueError(
                f"Unsupported data type: {self.data_type}. Supported types are 'visium' and 'xenium'."
            )
        hex_coordinates_dict = [
            {barcode: hex_coord}
            for hex_coord, barcode in zip(hex_coords, barcodes)
        ]
        coordinates_dict = {
            barcode: coord for barcode, coord in zip(barcodes, coords)
        }

        return hex_coordinates_dict, coordinates_dict

    def get_obsm(self, key, barcode, col=None, dtype=float):
        """
        Return 'adata.obsm[key].loc[barcode, col]' where col defaults to the
        first column of the dataframe.
        """
        if col is None:
            col = self.anndata.obsm[key].columns[0]  # Default to first column
        return dtype(self.anndata.obsm[key].loc[barcode, col])

    def get_X(self, barcode, gene=None, dtype=float):
        """
        Return 'adata[barcode, gene].X.toarray()' where gene defaults to the
        first gene in var.index. If the result is a single number, return as
        dtype.
        """
        if gene is None:
            gene = self.anndata.var.index[0]

        X = self.anndata[barcode, gene].X

        # Handle both sparse and dense matrices
        if sparse.issparse(X):
            expressions = X.toarray()
        else:
            expressions = np.asarray(X)

        #original_gene = self.var_upper_to_name.get(gene.upper())
        #expressions = self.anndata[barcode, original_gene].X.toarray()
        if expressions.size == 1:
            return dtype(expressions.flatten()[0])


    def to_geojson(self):
        hexagons = {"type": "FeatureCollection", "features": []}

        gene_expression_of_interest = None
        for score in genewise_scores:
            if score in self.anndata.uns:
                gene_expression_of_interest = self.anndata.uns[score].index[0]
                break  # Moran's I has prio over Geary's C

        for barcode, coords in zip(self.anndata.obs.index, self.coordinates):
            if (
                "in_tissue" in self.obs[barcode]
                and self.obs[barcode]["in_tissue"] == 0
            ):
                continue

            property_dict = {
                "barcode": barcode,
            }

            leiden_cluster = (
                int(self.anndata.obs.get("leiden", {}).get(barcode, -1))
                if "leiden" in self.anndata.obs.columns
                else None
            )

            if leiden_cluster is not None:
                property_dict["leiden"] = leiden_cluster

                if "leiden_centrality_scores" in self.anndata.uns:
                    property_dict["leiden_centrality"] = (
                        self.anndata.uns["leiden_centrality_scores"]
                        .iloc[leiden_cluster]
                        .to_dict()
                    )

                if "leiden_co_occurrence" in self.anndata.uns:
                    property_dict["leiden_co_occurrence"] = self.anndata.uns[
                        "leiden_co_occurrence"
                    ]["occ"][leiden_cluster].tolist()

                if "leiden_nhood_enrichment" in self.anndata.uns:
                    property_dict["leiden_nhood_enrichment"] = (
                        self.anndata.uns["leiden_nhood_enrichment"]["zscore"][
                            leiden_cluster
                        ].tolist()
                    )

            property_dict["centroid"] = (
                self.centers[barcode] if barcode in self.centers else None
            )


            for score in genie3_score_names + sponge_score_names:
                if score in self.anndata.obsm:
                    first_col = self.anndata.obsm[score].columns[0]
                    property_dict['regulatory_scores'] = self.anndata.obsm[
                        score
                    ].loc[barcode, first_col]
                    break


            score_mappings = {
                "ligand_receptor_relationships": (
                    "ligand_receptor_cosine_similarity",
                    "ligand_receptor_global_scores",
                ),
                "cell_comp_tf_activity_similarity": (
                    "cell_comp_tf_activity_cosine_similarity",
                    "cell_comp_tf_activity_global_scores",
                ),
                "tf_activity": ("tf_activity_score_ulm", None),
                "pathway_activity": ("pathway_activity_score_mlm", None),
            }

            for name, (obsm_key, uns_key) in score_mappings.items():
                if obsm_key in self.anndata.obsm:
                    col = (
                        self.anndata.uns[uns_key].index[0] if uns_key else None
                    )
                    property_dict[name] = self.get_obsm(obsm_key, barcode, col)

            for score in genewise_scores:
                if gene_expression_of_interest is not None:
                    property_dict["gene_expression"] = self.get_X(
                        barcode, gene_expression_of_interest
                    )
                    break


            # Add motif group scores if chromvar_activity is available
            if self.motif_groups and "chromvar_activity" in self.anndata.obsm:
                chromvar_df = self.anndata.obsm["chromvar_activity"]

                # case-insensitive mapping motif_name_upper -> actual column name
                motif_upper_to_col = {
                    m.upper(): m for m in chromvar_df.columns
                }

                # total sum of all motifs
                all_cols = [
                    motif_upper_to_col.get(motif_name.upper())
                    for motif_name in self.motif_names
                ]
                all_cols = [c for c in all_cols if c is not None]

                total_sum = float(chromvar_df.loc[barcode, all_cols].sum()) if all_cols else 0.0
                property_dict["chromvar_total_sum"] = total_sum

                for group_name, motif_list in self.motif_groups.items():
                    cols_for_group = []
                    for motif_name in motif_list:
                        col = motif_upper_to_col.get(motif_name.upper())
                        if col is not None:
                            cols_for_group.append(col)

                    if not cols_for_group:
                        # none of the motifs in this group found, skip
                        continue

                    group_score = float(
                        chromvar_df.loc[barcode, cols_for_group].sum()
                    )
                    property_dict[f"motif_{group_name}"] = group_score


            # Add additional properties from obs
            for key, value in self.obs[barcode].items():
                if value is None or value == "":
                    continue
                # Check for NaN if value is a float
                if isinstance(value, float) and np.isnan(value):
                    value = None
                property_dict[key] = value

            feature_dict = {
                "type": "Feature",
                "geometry": {
                    "type": self.geometry_type,
                    "coordinates": [coords[barcode]],
                },
                "properties": property_dict,
            }

            hexagons["features"].append(feature_dict)

        return hexagons


def load_adata(path: str) -> sc.AnnData:
    # NOTE (Xenium):
    # This script expects CELL-LEVEL AnnData with broadcasted grid scores
    # (e.g. xenium_cells_with_grid_scores.h5ad).
    # Passing a grid-level AnnData (st_grid.h5ad) will make Xenium look like Visium.

    spatial_data = sc.read_h5ad(path)

    # Reconstruct liana columns for placeholder fields in the geojson
    reconstruct_obsm_cols = {
        "ligand_receptor_cosine_similarity": "ligand_receptor",
        "cell_comp_tf_activity_cosine_similarity": "cell_comp_tf_activity",
    }
    for obsm_key, col_names_key in reconstruct_obsm_cols.items():
        if (
            obsm_key in spatial_data.obsm
            and "liana_columns" in spatial_data.uns
            and col_names_key in spatial_data.uns["liana_columns"]
        ):
            spatial_data.obsm[obsm_key] = pd.DataFrame(
                spatial_data.obsm[obsm_key],
                columns=spatial_data.uns["liana_columns"][col_names_key],
                index=spatial_data.obs_names,
            )
    if (
        "chromvar_spot_scores" in spatial_data.obsm
        and "chromvar_motifs" in spatial_data.uns
    ):
        motif_names = list(spatial_data.uns["chromvar_motifs"])
        spatial_data.obsm["chromvar_activity"] = pd.DataFrame(
            spatial_data.obsm["chromvar_spot_scores"],
            index=spatial_data.obs_names,
            columns=motif_names,
        )

    return spatial_data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Convert Visium data to GeoJSON format."
    )

    parser.add_argument(
        "--adata",
        "-a",
        type=str,
        required=True,
        help="Path to the input .h5ad file.",
    )
    parser.add_argument(
        "--radius", "-r", type=int, default=50, help="Radius of the hexagons."
    )
    parser.add_argument(
        "--scale",
        "-sc",
        type=float,
        default=0.5,
        help="Scale factor for the hexagons.",
    )
    parser.add_argument(
        "--data_type",
        "-dt",
        type=str,
        default="visium",
        choices=["visium", "xenium"],
        help="Type of spatial data (visium or xenium).",
    )
    parser.add_argument(
        "--outpath",
        "-o",
        type=str,
        default="./hexagons.geojson",
        help="Output path for the GeoJSON file.",
    )

    parser.add_argument(
        "--motif_groups",
        type=str,
        default=None,#"../backend/data/motif_groups.json",
        # make this a dict
            #     motif_groups = {
            #       "GATA_like": ["m1", "m4"],
            #       "RUNX": ["m3"],
            #       "MYC_cluster": ["m10", "m20"],
            #       }
        help="JSON dict: {group_name: [motif1, motif2, ...]} for chromVAR motif groups.",
    )


    args = parser.parse_args()
    spatial_data = load_adata(args.adata)

    # Sort global liana scores by cosine similarity std
    # We do this so that the tables appear sorted on the website
    global_scores_sort_keys = {
        "ligand_receptor_global_scores": "cosine_similarity_std",
        "cell_comp_tf_activity_global_scores": "cosine_similarity_std",
        "moranI": "I",
        "gearyC": "C",
        "chromvar_moranI": "I",
        "chromvar_gearyC": "C",
    }
    for global_score, sort_key in global_scores_sort_keys.items():
        if global_score in spatial_data.uns:
            spatial_data.uns[global_score] = spatial_data.uns[
                global_score
            ].sort_values(sort_key, ascending=False)


    if args.motif_groups is not None:
        motif_groups = json.loads(args.motif_groups)
    else:
        motif_groups = {"GATA_like": [
                            "MA0076.3",
                            "MA0079.5"
                        ],
                        "RUNX_like": [
                            "MA0003.5"
                        ],
                        "MYC_cluster": [
                            "MA0140.3",
                            "MA0141.4"
                        ],
                        "NFkB_like": [
                            "MA0105.4",
                            "MA0107.1"
                        ]
                        }


        #motif_groups = json.loads(args.motif_groups)

    hexagons = Hexagons(
        spatial_data,
        radius=args.radius,
        scale=args.scale,
        data_type=args.data_type,
        motif_groups=motif_groups,
    )

    print("uns keys:", sorted([k for k in spatial_data.uns.keys()
                          if "moran" in k.lower() or "geary" in k.lower()]))

    geojson_data = hexagons.to_geojson()

    # Add meta information like ligand receptor pair names for api fetching
    meta_dict = {}

    if any(gs in spatial_data.uns for gs in global_scores_sort_keys):
        for global_score in global_scores_sort_keys:
            if global_score in spatial_data.uns:
                meta_dict[global_score] = spatial_data.uns[global_score].to_dict()

    meta_dict['global_regulatory_scores_genie3'] = {score: spatial_data.obsm[score].mean().to_dict() for score in genie3_score_names + sponge_score_names if score in spatial_data.obsm and score.endswith('_genie3')}
    meta_dict['global_regulatory_scores_sponge'] = {score: spatial_data.obsm[score].mean().to_dict() for score in genie3_score_names + sponge_score_names if score in spatial_data.obsm and score.endswith('_sponge')}


    # The names in the tuple are options; all should have the same column names
    # but we don't want to rely on one obsm key being there
    colname_mapping = {
        "nmf_factors": ("ligand_receptor_NMF_factors",),
        "tf_names": ("tf_activity_score_ulm", "pathway_activity_padj_ulm"),
        "pathway_names": (
            "pathway_activity_score_mlm",
            "pathway_activity_padj_mlm",
        ),
        # "genie3_geneset_names": genie3_score_names,
        # "sponge_geneset_names": sponge_score_names,
    }
    for meta_key, obsm_key_options in colname_mapping.items():
        for obsm_key in obsm_key_options:
            if obsm_key in spatial_data.obsm:
                meta_dict[meta_key] = spatial_data.obsm[
                    obsm_key
                ].columns.tolist()
                break

    for geneset_key in ["genie_genesets", "sponge_genesets"]:
        if geneset_key in spatial_data.uns:
            genesets = {}
            for k, v in spatial_data.uns[geneset_key].items():
                genesets[k] = v.tolist()
            meta_dict[geneset_key] = genesets

    available_grn_score_names = [
        grn_score_name
        for grn_score_name in genie3_score_names + sponge_score_names
        if grn_score_name in spatial_data.obsm
    ]

    if len(available_grn_score_names) > 0:
        meta_dict["grn_score_names"] = available_grn_score_names

    if (
        "leiden_co_occurrence" in spatial_data.uns
        and "interval" in spatial_data.uns["leiden_co_occurrence"]
    ):
        meta_dict["interval"] = spatial_data.uns["leiden_co_occurrence"][
            "interval"
        ].tolist()

    meta_dict["data_type"] = args.data_type
    # Add differential motif activity top-motif tables
    if "diff_motif_activity_top_motifs" in spatial_data.uns:
        meta_dict["diff_motif_activity_top_motifs"] = {}

        for comparison, df in spatial_data.uns["diff_motif_activity_top_motifs"].items():
            # If stored as DataFrame, convert to table-friendly dict
            if isinstance(df, pd.DataFrame):
                meta_dict["diff_motif_activity_top_motifs"][comparison] = df.to_dict()
            else:
                # fallback in case it already came in as a plain dict-like object
                meta_dict["diff_motif_activity_top_motifs"][comparison] = pd.DataFrame(df).to_dict()

    geojson_data["meta"] = meta_dict

    os.makedirs(os.path.dirname(args.outpath), exist_ok=True)

    with open(args.outpath, "w+") as f:
        json.dump(geojson_data, f, indent=4, ignore_nan=True)

    print("GeoJSON file created successfully.")
