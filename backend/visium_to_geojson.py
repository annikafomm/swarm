import pandas as pd
import numpy as np
import scanpy as sc
import geopandas as gpd
import json
import argparse
import os


class Hexagons:
    """
    Class to convert Visium or Transformed Xenium spatial data to GeoJSON format with hexagonal geometry.
    This class handles the parsing of AnnData objects, Spatial scores, Gene regulatory Scores and spatial coordinates
    to create a GeoJSON representation of the spatial data with annotated properties.
    """


    def __init__(self, anndata, aucell_genie3_path, aucell_sponge_path, radius=5, scale=0.1, data_type="visium"):

        self.anndata = anndata
        if aucell_genie3_path and aucell_sponge_path:
            if not os.path.exists(aucell_genie3_path):
                raise FileNotFoundError(f"Aucell Genie3 CSV file not found: {aucell_genie3_path}")
            if not os.path.exists(aucell_sponge_path):
                raise FileNotFoundError(f"Aucell Sponge CSV file not found: {aucell_sponge_path}")
            self.parse_aucell_csvs(aucell_genie3_path, aucell_sponge_path)
        self.radius = radius
        self.type = "FeatureCollection"
        self.geometry_type = "Polygon"
        self.scale = scale
        self.data_type = data_type
        self.coordinates, self.centers = self.parse_coordinates()
        self.obs = self.parse_obs()



    def parse_obs(self):
        obs = self.anndata.obs
        properties = obs.copy()
        # Convert to dict with id as key and properties as values
        properties_dict = properties.to_dict(orient='index')
        return properties_dict


    """
    def parse_aucell_gene_sets(self):
        uns = self.anndata.uns
        gene_sets = uns.get("genesets", None)
        if gene_sets is not None:
            return gene_sets
        else:
            print("No Aucell gene sets found in the AnnData object.")
            return {}
    Methods are parse from csv files, so this is no longer used and will be removed in future commits.
    """

    def hexagon_points(self,x, y, radius):
        # Scale centers
        x = x * self.scale
        y = y * self.scale
        return [
            (x + radius * np.cos(np.pi / 3 * i), y + radius * np.sin(np.pi / 3 * i))
            for i in range(6)
        ] + [(x + radius * np.cos(0), y + radius * np.sin(0))]

    def parse_coordinates(self):
        if self.data_type == "visium":
            anndata_spatial_coordinates = self.anndata.obsm["spatial"].copy()
            hex_coords = []
            coords = []
            for coord_tuple in anndata_spatial_coordinates:
                x, y = coord_tuple
                hexagon = self.hexagon_points(x=x, y=y, radius=self.radius)
                hex_coords.append(hexagon)
                coords.append([int(x) * self.scale, int(y) * self.scale])
            barcodes = self.anndata.obs.index
            # Convert to list of tuples (x, y)
            hex_coordinates_dict = [{barcode: hex_coord} for hex_coord, barcode in zip(hex_coords, barcodes)]
            coordinates_dict = {barcode: coord for barcode, coord in zip(barcodes, coords)}
            return hex_coordinates_dict, coordinates_dict
        if self.data_type == "xenium":
            print("Xenium data type is not implemented yet.")
            return [], {}
        else:
            raise ValueError(f"Unsupported data type: {self.data_type}. Supported types are 'visium' and 'xenium'.")

    def parse_aucell_csvs(self, aucell_genie3_path, aucell_sponge_path):
        print("Parsing Aucell CSV files...")
        aucell_scores_genie3 = pd.read_csv(aucell_genie3_path, index_col=0)
        aucell_scores_sponge = pd.read_csv(aucell_sponge_path, index_col=0)
        self.aucell_scores_genie3 = aucell_scores_genie3
        self.aucell_scores_sponge = aucell_scores_sponge
        print("Aucell CSV files parsed successfully.")


    def to_geojson(self):
        hexagons = {
                "type": "FeatureCollection",
                "features": []
            }
        for barcode, coords in zip(self.anndata.obs.index, self.coordinates):

            if self.obs[barcode]["in_tissue"] == False:
                continue

            leiden_cluster = int(self.anndata.obs.get("leiden", {}).get(barcode, None))

            leiden_centrality = self.anndata.uns.get("leiden_centrality_scores", {}).iloc[leiden_cluster].to_dict()
            leiden_co_occurrence = self.anndata.uns.get("leiden_co_occurrence")["occ"][leiden_cluster]
            leiden_nhood_zscores = self.anndata.uns.get("leiden_nhood_enrichment",{}).get("zscore")[leiden_cluster]


            # Print all types
            print(type(leiden_centrality))
            print(type(leiden_co_occurrence))
            print(type(leiden_nhood_zscores))

            feature_dict = {"type": "Feature",
                "geometry": {
                    "type": self.geometry_type,
                    "coordinates": [coords[barcode]]
                },
                "properties": {
                    "barcode": barcode,
                    "aucell_genie3": self.aucell_scores_genie3.get(barcode, {}).to_dict() if barcode in self.aucell_scores_genie3 else {},
                    "aucell_sponge": self.aucell_scores_sponge.get(barcode, {}).to_dict() if barcode in self.aucell_scores_sponge else {},
                    "leiden": leiden_cluster,
                    "leiden_centrality": leiden_centrality,
                    "leiden_co_occurrence": leiden_co_occurrence.tolist(),
                    "leiden_nhood_enrichment": leiden_nhood_zscores.tolist(),
                    "centroid": self.centers[barcode] if barcode in self.centers else None
                    }
                }

            # Add additional properties from obs
            for key, value in self.obs[barcode].items():

                if value is None or value == "":
                    continue
            # Check for NaN if value is a float
                if isinstance(value, float) and np.isnan(value):
                    continue

                feature_dict["properties"][key] = value

            # Add leiden properties from uns--



            hexagons["features"].append(feature_dict)

        return hexagons



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Visium data to GeoJSON format.")

    parser.add_argument("--adata","-a", type=str, required=True, help="Path to the input .h5ad file.")
    parser.add_argument("--aucell_genie3_path","-g", type=str, required=True, help="Path to the Aucell Genie3 CSV file.")
    parser.add_argument("--aucell_sponge_path", "-s", type=str, required=True, help="Path to the Aucell Sponge CSV file.")
    parser.add_argument("--radius","-r", type=int, default=50, help="Radius of the hexagons.")
    parser.add_argument("--scale", "-sc", type=float, default=0.5, help="Scale factor for the hexagons.")
    parser.add_argument("--data_type", "-dt", type=str, default="visium", choices=["visium", "xenium"], help="Type of spatial data (visium or xenium).")
    parser.add_argument("--outpath", "-o", type=str, default="./hexagons.geojson", help="Output path for the GeoJSON file.")
    args = parser.parse_args()

    spatial_data = sc.read_h5ad(args.adata)

    hexagons = Hexagons(spatial_data, radius=args.radius, scale=0.4, data_type=args.data_type, aucell_genie3_path=args.aucell_genie3_path, aucell_sponge_path=args.aucell_sponge_path)

    geojson_data = hexagons.to_geojson()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert Visium data to GeoJSON format.")

    parser.add_argument("--adata","-a", type=str, required=True, help="Path to the input .h5ad file.")
    parser.add_argument("--aucell_genie3_path","-g", type=str, required=False, help="Path to the Aucell Genie3 CSV file.")
    parser.add_argument("--aucell_sponge_path", "-s", type=str, required=False, help="Path to the Aucell Sponge CSV file.")
    parser.add_argument("--radius","-r", type=int, default=50, help="Radius of the hexagons.")
    parser.add_argument("--scale", "-sc", type=float, default=0.5, help="Scale factor for the hexagons.")
    parser.add_argument("--data_type", "-dt", type=str, default="visium", choices=["visium", "xenium"], help="Type of spatial data (visium or xenium).")
    parser.add_argument("--outpath", "-o", type=str, default="./hexagons.geojson", help="Output path for the GeoJSON file.") # Files should be named after their dataset names so you can select them in the view
    args = parser.parse_args()

    spatial_data = sc.read_h5ad(args.adata)

    if not args.aucell_genie3_path or not args.aucell_sponge_path:
        hexagons = Hexagons(spatial_data, radius=args.radius, scale=args.scale, data_type=args.data_type)
    else:
        hexagons = Hexagons(spatial_data, radius=args.radius, scale=args.scale, data_type=args.data_type, aucell_genie3_path=args.aucell_genie3_path, aucell_sponge_path=args.aucell_sponge_path)

    geojson_data = hexagons.to_geojson()

    os.makedirs(os.path.dirname(args.outpath), exist_ok=True)

    with open(args.outpath, "w+") as f:
        json.dump(geojson_data, f, indent=4)

    print("GeoJSON file created successfully.")


