import { Component, EventEmitter, Input, OnChanges, Output, SimpleChanges } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { MatIconModule } from '@angular/material/icon';
import { MatButtonModule } from '@angular/material/button';
import { MatTooltipModule } from '@angular/material/tooltip';
import { CellFeature } from '../hexagon-view/cell-feature.types';
import { InfoService } from '../info.service';
import { TranslationService } from '../translation.service';

export interface AttributeStats {
  key: string;
  label: string;
  type: 'numeric' | 'categorical';
  uniqueCount: number;
  count: number;
  min: number | null;
  max: number | null;
  mean: number | null;
  variance: number | null;
  moranI: number | null;
  getisOrdG: number | null;
  getisOrdZ: number | null;
  topCategory?: string;
  topCategoryCount?: number;
  description?: string | null;
}

@Component({
  selector: 'app-further-attributes-panel',
  standalone: true,
  imports: [CommonModule, FormsModule, MatIconModule, MatButtonModule, MatTooltipModule],
  templateUrl: './further-attributes-panel.component.html',
  styleUrls: ['./further-attributes-panel.component.scss'],
})
export class FurtherAttributesPanelComponent implements OnChanges {
  @Input() features: CellFeature[] = [];
  @Input() selectedAttribute: string | null = null;
  @Input() isCompare: boolean = false;
  @Output() attributeSelected = new EventEmitter<string>();

  attributes: AttributeStats[] = [];
  filteredAttributes: AttributeStats[] = [];
  searchQuery: string = '';
  typeFilter: 'all' | 'numeric' | 'categorical' = 'all';

  sortColumn: keyof AttributeStats = 'label';
  sortAsc: boolean = true;

  // Pagination
  pageSize: number = 10;
  currentPage: number = 1;

  // Spatial neighbor graph cache
  private neighborIndices: number[][] = [];
  private lastFeaturesRef: CellFeature[] | null = null;

  // Known properties mapped to other dedicated tabs
  private readonly reservedTabKeys = new Set<string>([
    'leiden',
    'leiden_centrality',
    'leiden_co_occurrence',
    'leiden_nhood_enrichment',
    'degree_centrality',
    'average_clustering',
    'closeness_centrality',
    'co_occurrence',
    'regulatory_scores',
    'aucell_genie3',
    'aucell_sponge',
    'ligand_receptor_relationships',
    'cell_comp_tf_activity_similarity',
    'tf_activity',
    'pathway_activity',
    'gene_expression',
    'chromvar_total_sum',
    'chromvar_moranI',
    'chromvar_gearyC',
    'diff_motif_activity_top_motifs',
    'color',
  ]);

  constructor(
    public infoService: InfoService,
    private translationService: TranslationService,
  ) {}

  ngOnChanges(changes: SimpleChanges): void {
    if (changes['features']) {
      this.computeAttributes();
    }
  }

  /**
   * Builds the spatial nearest-neighbor graph ($k=6$) and computes stats for all eligible attributes.
   */
  private computeAttributes(): void {
    if (!this.features || this.features.length === 0) {
      this.attributes = [];
      this.filteredAttributes = [];
      this.neighborIndices = [];
      return;
    }

    const n = this.features.length;

    // 1. Build or reuse spatial neighbor graph if features changed
    if (this.features !== this.lastFeaturesRef) {
      this.lastFeaturesRef = this.features;
      this.neighborIndices = this.buildSpatialNeighbors(this.features, 6);
    }

    // 2. Discover all unique keys across all features
    const allKeys = new Set<string>();
    for (let i = 0; i < n; i++) {
      const props = this.features[i]?.properties;
      if (props) {
        for (const k of Object.keys(props)) {
          allKeys.add(k);
        }
      }
    }

    const statsList: AttributeStats[] = [];

    for (const key of allKeys) {
      const lowerKey = key.toLowerCase();

      // Skip keys that belong to dedicated tabs or internal state
      if (this.reservedTabKeys.has(lowerKey) || this.reservedTabKeys.has(key)) {
        continue;
      }
      if (lowerKey.endsWith('_genie3') || lowerKey.endsWith('_sponge')) {
        continue;
      }

      // Collect values across all features
      const values: unknown[] = [];
      for (let i = 0; i < n; i++) {
        const val = this.features[i]?.properties?.[key];
        if (val !== undefined && val !== null && val !== '' && String(val).toLowerCase() !== 'nan') {
          values.push(val);
        }
      }

      const validCount = values.length;
      if (validCount === 0) continue;

      // Count unique values
      const uniqueSet = new Set<string>();
      for (let i = 0; i < validCount; i++) {
        const v = values[i];
        uniqueSet.add(typeof v === 'object' ? JSON.stringify(v) : String(v));
      }
      const uniqueCount = uniqueSet.size;

      // Heuristic 1: Constant values (all cells share 1 identical value e.g. patient ID, treatment, assay) -> skip
      if (uniqueCount <= 1) {
        continue;
      }

      // Heuristic 2: Unique cell identifier or coordinate (every cell has unique distinct ID e.g. barcode, cell_id) -> skip
      if (
        (uniqueCount === validCount && (typeof values[0] === 'string' || validCount > 20)) ||
        (validCount > 50 && uniqueCount >= validCount * 0.98 && typeof values[0] === 'string') ||
        ['barcode', 'centroid', 'observation_joinid', 'cell_id', 'id', 'guid'].includes(lowerKey)
      ) {
        continue;
      }

      // Determine if numeric
      const isAllNumeric = values.every((v) => {
        if (typeof v === 'number') return Number.isFinite(v);
        if (typeof v === 'string') return v.trim() !== '' && !Number.isNaN(Number(v));
        return false;
      });

      const label = this.formatLabel(key);
      const desc = this.infoService.getFieldInfo(key);

      if (isAllNumeric) {
        const numValues = new Float64Array(n);
        let filledCount = 0;
        let sum = 0;
        let min = Infinity;
        let max = -Infinity;

        for (let i = 0; i < n; i++) {
          const raw = this.features[i]?.properties?.[key];
          if (raw !== undefined && raw !== null && raw !== '' && !Number.isNaN(Number(raw))) {
            const num = Number(raw);
            numValues[i] = num;
            sum += num;
            if (num < min) min = num;
            if (num > max) max = num;
            filledCount++;
          } else {
            numValues[i] = NaN;
          }
        }

        if (filledCount < 2) continue;

        const mean = sum / filledCount;

        // Compute variance
        let sumSqDiff = 0;
        for (let i = 0; i < n; i++) {
          const val = numValues[i];
          if (!Number.isNaN(val)) {
            const diff = val - mean;
            sumSqDiff += diff * diff;
          }
        }

        const variance = filledCount > 1 ? sumSqDiff / (filledCount - 1) : 0;
        if (variance <= 1e-12) {
          // No numeric variance -> skip
          continue;
        }

        // Compute Moran's I and Getis-Ord G using spatial neighbor graph
        const { moranI, getisOrdG, getisOrdZ } = this.computeSpatialMetrics(numValues, mean, sumSqDiff);

        statsList.push({
          key,
          label,
          type: 'numeric',
          uniqueCount,
          count: filledCount,
          min,
          max,
          mean,
          variance,
          moranI,
          getisOrdG,
          getisOrdZ,
          description: desc,
        });
      } else {
        // Categorical attribute (e.g. cell type annotation, region, cluster call)
        // High-cardinality string IDs that are essentially random -> limit to reasonable categories
        if (uniqueCount > Math.min(60, Math.floor(n * 0.5))) {
          continue;
        }

        // Find most common category
        const counts = new Map<string, number>();
        for (let i = 0; i < validCount; i++) {
          const str = String(values[i]);
          counts.set(str, (counts.get(str) || 0) + 1);
        }

        let topCat = '';
        let topCount = 0;
        for (const [cat, cnt] of counts.entries()) {
          if (cnt > topCount) {
            topCount = cnt;
            topCat = cat;
          }
        }

        statsList.push({
          key,
          label,
          type: 'categorical',
          uniqueCount,
          count: validCount,
          min: null,
          max: null,
          mean: null,
          variance: null,
          moranI: null,
          getisOrdG: null,
          getisOrdZ: null,
          topCategory: topCat,
          topCategoryCount: topCount,
          description: desc,
        });
      }
    }

    this.attributes = statsList;
    this.applyFilters();
  }

  /**
   * Fast spatial nearest neighbors calculation ($k=6$) using grid partitioning or pairwise search.
   */
  private buildSpatialNeighbors(features: CellFeature[], k: number = 6): number[][] {
    const n = features.length;
    const coords: [number, number][] = new Array(n);

    for (let i = 0; i < n; i++) {
      coords[i] = this.extractCentroid(features[i]);
    }

    const neighbors: number[][] = new Array(n);

    // For typical Visium/Xenium datasets (N up to ~10,000), spatial grid binning gives O(N) neighbor lookup
    const minX = Math.min(...coords.map((c) => c[0]));
    const maxX = Math.max(...coords.map((c) => c[0]));
    const minY = Math.min(...coords.map((c) => c[1]));
    const maxY = Math.max(...coords.map((c) => c[1]));

    const gridSize = Math.max(1, Math.sqrt(n) / 2);
    const cellW = (maxX - minX || 1) / gridSize;
    const cellH = (maxY - minY || 1) / gridSize;

    const grid = new Map<string, number[]>();
    for (let i = 0; i < n; i++) {
      const gx = Math.floor((coords[i][0] - minX) / cellW);
      const gy = Math.floor((coords[i][1] - minY) / cellH);
      const gkey = `${gx},${gy}`;
      const bucket = grid.get(gkey);
      if (bucket) {
        bucket.push(i);
      } else {
        grid.set(gkey, [i]);
      }
    }

    for (let i = 0; i < n; i++) {
      const [xi, yi] = coords[i];
      const gx = Math.floor((xi - minX) / cellW);
      const gy = Math.floor((yi - minY) / cellH);

      // Collect candidates from adjacent grid cells
      const candidates: number[] = [];
      for (let dx = -2; dx <= 2; dx++) {
        for (let dy = -2; dy <= 2; dy++) {
          const b = grid.get(`${gx + dx},${gy + dy}`);
          if (b) {
            for (let idx = 0; idx < b.length; idx++) {
              if (b[idx] !== i) {
                candidates.push(b[idx]);
              }
            }
          }
        }
      }

      // If grid candidate count is too small (e.g. isolated boundary), fallback to all
      const searchPool = candidates.length >= k ? candidates : Array.from({ length: n }, (_, idx) => idx).filter((idx) => idx !== i);

      // Compute squared Euclidean distances to candidates
      const dists: { index: number; distSq: number }[] = new Array(searchPool.length);
      for (let j = 0; j < searchPool.length; j++) {
        const cidx = searchPool[j];
        const dx = xi - coords[cidx][0];
        const dy = yi - coords[cidx][1];
        dists[j] = { index: cidx, distSq: dx * dx + dy * dy };
      }

      dists.sort((a, b) => a.distSq - b.distSq);
      neighbors[i] = dists.slice(0, k).map((d) => d.index);
    }

    return neighbors;
  }

  /**
   * Extracts or computes spatial centroid [x, y] from any CellFeature geometry or properties.
   */
  private extractCentroid(f: CellFeature): [number, number] {
    const props: any = f.properties || {};
    if (Array.isArray(props['centroid']) && props['centroid'].length >= 2) {
      return [Number(props['centroid'][0]), Number(props['centroid'][1])];
    }
    if (props['x'] !== undefined && props['y'] !== undefined) {
      return [Number(props['x']), Number(props['y'])];
    }
    if (props['spatial_x'] !== undefined && props['spatial_y'] !== undefined) {
      return [Number(props['spatial_x']), Number(props['spatial_y'])];
    }
    if (props['array_col'] !== undefined && props['array_row'] !== undefined) {
      return [Number(props['array_col']) * 10, Number(props['array_row']) * 10];
    }

    // Geometry coordinates
    const geom = f.geometry as any;
    if (geom) {
      if (geom.type === 'Point' && Array.isArray(geom.coordinates)) {
        return [Number(geom.coordinates[0]), Number(geom.coordinates[1])];
      }
      if (Array.isArray(geom.coordinates?.[0])) {
        const ring = geom.coordinates[0];
        let sx = 0;
        let sy = 0;
        const len = ring.length;
        for (let i = 0; i < len; i++) {
          sx += Number(ring[i][0]);
          sy += Number(ring[i][1]);
        }
        return [sx / len, sy / len];
      }
    }

    return [0, 0];
  }

  /**
   * Computes Moran's I and Getis-Ord G statistics from values and neighbor adjacency graph.
   */
  private computeSpatialMetrics(
    values: Float64Array,
    mean: number,
    sumSqDiff: number,
  ): { moranI: number | null; getisOrdG: number | null; getisOrdZ: number | null } {
    const n = values.length;
    if (n < 4 || sumSqDiff <= 1e-12) {
      return { moranI: null, getisOrdG: null, getisOrdZ: null };
    }

    let moranNumerator = 0;
    let getisNumerator = 0;
    let sumX = 0;
    let sumX2 = 0;
    let totalPairs = 0;

    for (let i = 0; i < n; i++) {
      const xi = values[i];
      if (Number.isNaN(xi)) continue;

      sumX += xi;
      sumX2 += xi * xi;

      const nbrs = this.neighborIndices[i] || [];
      const k = nbrs.length;
      if (k === 0) continue;

      const diffI = xi - mean;
      let nbrDiffSum = 0;
      let nbrValSum = 0;
      let validNbrs = 0;

      for (let j = 0; j < k; j++) {
        const xj = values[nbrs[j]];
        if (!Number.isNaN(xj)) {
          nbrDiffSum += xj - mean;
          nbrValSum += xj;
          validNbrs++;
        }
      }

      if (validNbrs > 0) {
        moranNumerator += diffI * (nbrDiffSum / validNbrs);
        getisNumerator += xi * nbrValSum;
        totalPairs += validNbrs;
      }
    }

    const moranI = sumSqDiff > 0 ? moranNumerator / sumSqDiff : null;

    // Getis-Ord General G
    const denomG = sumX * sumX - sumX2;
    let getisOrdG: number | null = null;
    let getisOrdZ: number | null = null;

    if (denomG > 0 && totalPairs > 0) {
      getisOrdG = getisNumerator / denomG;

      // Theoretical expectation and variance of G under null hypothesis of spatial randomness
      const expectedG = totalPairs / (n * (n - 1));
      const s1 = totalPairs * 2; // Symmetric binary edges
      const s2 = 4 * totalPairs;
      const b0 = (n * n - 3 * n + 3) * s1 - n * s2 + 3 * totalPairs * totalPairs;
      const b1 = -((n * n - n) * s1 - 2 * n * s2 + 6 * totalPairs * totalPairs);
      const b2 = -(2 * n * s1 - (n + 3) * s2 + 6 * totalPairs * totalPairs);
      const b3 = 4 * (n - 1) * s1 - 2 * (n + 1) * s2 + 8 * totalPairs * totalPairs;
      const b4 = s1 - s2 + (totalPairs * totalPairs);

      const m1 = sumX;
      const m2 = sumX2;
      const m3 = values.reduce((acc, v) => acc + (Number.isNaN(v) ? 0 : v * v * v), 0);
      const m4 = values.reduce((acc, v) => acc + (Number.isNaN(v) ? 0 : v * v * v * v), 0);

      const d = (m1 * m1 - m2) * (m1 * m1 - m2);
      const eg2Num = (b0 * (m2 * m2) + b1 * (m1 * m3) + b2 * (m1 * m1 * m2) + b3 * (m1 * m1 * m1 * m1) + b4 * m4);
      const eg2Denom = (n * (n - 1) * (n - 2) * (n - 3) * d) / (n * n);
      const varG = eg2Denom > 0 ? Math.max(0, (eg2Num / eg2Denom) - expectedG * expectedG) : 0;

      if (varG > 1e-12) {
        getisOrdZ = (getisOrdG - expectedG) / Math.sqrt(varG);
      }
    }

    return { moranI, getisOrdG, getisOrdZ };
  }

  private formatLabel(key: string): string {
    const translated = this.translationService.translateSync(key);
    if (translated !== key) return translated;
    return key.replace(/_/g, ' ').replace(/\b\w/g, (c) => c.toUpperCase());
  }

  applyFilters(): void {
    const query = this.searchQuery.trim().toLowerCase();

    this.filteredAttributes = this.attributes.filter((attr) => {
      // Type filter
      if (this.typeFilter !== 'all' && attr.type !== this.typeFilter) {
        return false;
      }

      // Search query
      if (query) {
        const matchesKey = attr.key.toLowerCase().includes(query);
        const matchesLabel = attr.label.toLowerCase().includes(query);
        const matchesTop = attr.topCategory?.toLowerCase().includes(query);
        return matchesKey || matchesLabel || !!matchesTop;
      }

      return true;
    });

    this.sortData();
    this.currentPage = 1;
  }

  setSort(col: keyof AttributeStats): void {
    if (this.sortColumn === col) {
      this.sortAsc = !this.sortAsc;
    } else {
      this.sortColumn = col;
      this.sortAsc = col === 'label' || col === 'key';
    }
    this.sortData();
  }

  private sortData(): void {
    const col = this.sortColumn;
    const asc = this.sortAsc ? 1 : -1;

    this.filteredAttributes.sort((a, b) => {
      const valA = a[col];
      const valB = b[col];

      if (valA === null || valA === undefined) return 1;
      if (valB === null || valB === undefined) return -1;

      if (typeof valA === 'number' && typeof valB === 'number') {
        return (valA - valB) * asc;
      }
      return String(valA).localeCompare(String(valB), undefined, { numeric: true }) * asc;
    });
  }

  selectAttribute(attr: AttributeStats): void {
    this.attributeSelected.emit(attr.key);
  }

  get paginatedAttributes(): AttributeStats[] {
    const start = (this.currentPage - 1) * this.pageSize;
    return this.filteredAttributes.slice(start, start + this.pageSize);
  }

  get totalPages(): number {
    return Math.max(1, Math.ceil(this.filteredAttributes.length / this.pageSize));
  }

  prevPage(): void {
    if (this.currentPage > 1) this.currentPage--;
  }

  nextPage(): void {
    if (this.currentPage < this.totalPages) this.currentPage++;
  }
}
