import pandas as pd
import numpy as np
import scanpy as sc
import squidpy as sq
import geopandas as gpd


class Hexagons:
    def __init__(self, anndata, radius=5, scale=0.1, data_type="visium"):

        self.anndata = anndata
        self.radius = radius
        self.type = "FeatureCollection"
        self.geometry_type = "Polygon"
        self.scale = scale
        self.data_type = data_type
        self.coordinates, self.centers = self.parse_coordinates()
        self.properties = self.parse_properties()
        
       
    
    def parse_properties(self):
        obs = self.anndata.obs
        properties = obs.copy()
        # Convert to dict with id as key and properties as values
        properties_dict = properties.to_dict(orient='index')
        return properties_dict
    
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


    
    def to_geojson(self):
        hexagons = {
                "type": "FeatureCollection",
                "features": []
            }
        for barcode, coords in zip(self.anndata.obs.index, self.coordinates):

            if self.properties[barcode]["in_tissue"] == False:
                continue

            hexagons["features"].append(
                {"type": "Feature",
                "geometry": {
                    "type": self.geometry_type,
                    "coordinates": [coords[barcode]]
                },
                "properties": {
                    "barcode": barcode,
                    "centroid": self.centers[barcode],
                    "cell_type": self.properties[barcode].get("cell_type", "unknown")
                    }
                }
            )
        
        return hexagons
    

class Counts:
    def __init__(self, anndata):
        self.anndata = anndata
        self.counts = self.parse_counts()

    def parse_counts(self):
        counts = self.anndata.X.toarray() if hasattr(self.anndata.X, "toarray") else self.anndata.X
        genes = self.anndata.var.index.tolist()
        barcodes = self.anndata.obs.index.tolist()
        counts_df = pd.DataFrame(counts, index=barcodes, columns=genes)
        return counts_df
        
    
if __name__ == "__main__":
    import json
    import argparse
    parser = argparse.ArgumentParser(description="Convert Visium data to GeoJSON format.")

    parser.add_argument("--adata", type=str, required=True, help="Path to the input .h5ad file.")
    parser.add_argument("--radius", type=int, default=5, help="Radius of the hexagons.")
    parser.add_argument("--scale", type=float, default=0.1, help="Scale factor for the hexagons.")
    parser.add_argument("--data_type", type=str, default="visium", choices=["visium", "xenium"], help="Type of spatial data (visium or xenium).")
    args = parser.parse_args()

    spatial_data = sc.read_h5ad(args.adata)

    hexagons = Hexagons(spatial_data, radius=args.radius, scale=0.4, data_type=args.data_type)
    counts = Counts(spatial_data)

    geojson_data = hexagons.to_geojson()

    with open("hexagons_GSM6592049_M2.geojson", "w") as f:
        json.dump(geojson_data, f, indent=4)
        
    print("GeoJSON file created successfully.")
        

        