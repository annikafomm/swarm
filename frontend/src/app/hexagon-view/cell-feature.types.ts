export interface CellGeometry {
  type: 'Polygon';
  coordinates: number[][][];
}

export interface CellProperties {
  barcode: string;
  centroid: [number, number] | [];
  cell_type: string;
  leiden: number;
  color: string;
  aucell_genie3: { [key: string]: number };
  aucell_sponge: { [key: string]: number };
  [key: string]: string | number | number[] | [] | undefined | { [key: string]: any };
}

export interface CellFeature {
  type: 'Feature';
  geometry: CellGeometry;
  properties: CellProperties;
}

export interface PropertyGroupItem {
  key: string;
  label: string;
  value: unknown;
  info: string | null;
}

export interface PropertyGroup {
  key: string;
  title: string;
  icon: string;
  description: string;
  items: PropertyGroupItem[];
}
