import { CellFeature } from './hexagon-view/cell-feature.types';

/**
 * Shared heuristics for deciding which dataset attributes are worth surfacing to the user (as
 * colorable/clusterable properties, in the Further Attributes table, etc). Extracted from
 * FurtherAttributesPanelComponent so every attribute-picking UI applies the same exclusions
 * instead of drifting apart.
 */

const KNOWN_NON_ATTRIBUTE_KEYS = new Set<string>([
  'barcode',
  'centroid',
  'observation_joinid',
  'cell_id',
  'id',
  'guid',
  'array_row',
  'array_col',
  'arrayrow',
  'arraycol',
  'array_x',
  'array_y',
  'array row',
  'array col',
  'spatial_x',
  'spatial_y',
  'spatial_row',
  'spatial_col',
  'x_original',
  'y_original',
  'original_x',
  'original_y',
  'orig_x',
  'orig_y',
  'x_orig',
  'y_orig',
  '_x',
  '_y',
  '__x',
  '__y',
  'x',
  'y',
  'pxl_col_in_fullres',
  'pxl_row_in_fullres',
  'pxl_col',
  'pxl_row',
  'pixel_x',
  'pixel_y',
  'imagecol',
  'imagerow',
  'image_col',
  'image_row',
  'image_x',
  'image_y',
  'x_coord',
  'y_coord',
  'coord_x',
  'coord_y',
  'coords_x',
  'coords_y',
  'center_x',
  'center_y',
  'centroid_x',
  'centroid_y',
  'x_centroid',
  'y_centroid',
  'grid_x',
  'grid_y',
  'grid_row',
  'grid_col',
  'spot_x',
  'spot_y',
  'spot_row',
  'spot_col',
]);

/**
 * Identifies keys that represent spatial coordinates, grid indices, pixel positions, or record
 * identifiers (barcode, cell_id, ...) — never meaningful to color, filter, or cluster by.
 */
export function isSpatialOrIdentifierKey(key: string): boolean {
  if (!key) return false;
  const lower = key.toLowerCase().trim();

  if (KNOWN_NON_ATTRIBUTE_KEYS.has(lower)) {
    return true;
  }

  if (/^_+[xy]$/i.test(lower) || /^[xy]$/i.test(lower)) {
    return true;
  }

  if (/^(?:spatial|array|grid|pxl|pixel|image|orig(?:inal)?|center|centroid|spot|coord|pos)_(?:x|y|row|col)$/i.test(lower)) {
    return true;
  }

  if (/^(?:x|y)_(?:orig(?:inal)?|coord|centroid|spatial|grid|pixel|pxl|pos|index)$/i.test(lower)) {
    return true;
  }

  return false;
}

/**
 * True when `key` has ≤1 unique value across `features` (constant, e.g. patient ID or assay), or
 * is a near-unique string identifier (e.g. barcode, cell_id) — both make it useless to color,
 * filter, or cluster by.
 */
export function isUninformativeAttribute(key: string, features: CellFeature[] | null | undefined): boolean {
  if (!features || features.length === 0) return false;

  const values: unknown[] = [];
  for (const f of features) {
    const val = f?.properties?.[key];
    if (val !== undefined && val !== null && val !== '' && String(val).toLowerCase() !== 'nan') {
      values.push(val);
    }
  }

  const validCount = values.length;
  if (validCount === 0) return true;

  const uniqueSet = new Set<string>();
  for (const v of values) {
    uniqueSet.add(typeof v === 'object' ? JSON.stringify(v) : String(v));
  }
  const uniqueCount = uniqueSet.size;

  if (uniqueCount <= 1) return true;

  const isAllNumeric = values.every((v) => {
    if (typeof v === 'number') return Number.isFinite(v);
    if (typeof v === 'string') return v.trim() !== '' && !Number.isNaN(Number(v));
    return false;
  });

  if (!isAllNumeric) {
    if (uniqueCount === validCount || (validCount > 50 && uniqueCount >= validCount * 0.98)) {
      return true;
    }
  }

  return false;
}
