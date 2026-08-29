import { Component, ElementRef, EventEmitter, Input, OnChanges, OnInit, Output, SimpleChanges, ViewChild } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { HttpClient } from '@angular/common/http';
import { MatButtonModule } from '@angular/material/button';
import { MatButtonToggleModule } from '@angular/material/button-toggle';
import { MatIconModule } from '@angular/material/icon';
import { MatTooltipModule } from '@angular/material/tooltip';
import { MatProgressSpinnerModule } from '@angular/material/progress-spinner';
import * as Plotly from 'plotly.js-dist-min';
import * as d3 from 'd3';
import { CellFeature } from '../hexagon-view/cell-feature.types';

export interface SpatialFeature {
  id: string;
  name: string;
  category: string;
  symbol?: string;
}

export interface PairScatterData {
  feature_a: SpatialFeature;
  feature_b: SpatialFeature;
  stats: {
    pearson_r: number;
    pearson_p: number;
    spearman_rho: number;
    r_squared: number;
    slope: number;
    intercept: number;
    bivariate_moran_i: number;
    is_dependent: boolean;
    dependency_reason: string;
    total_cells: number;
  };
  points: { cell_id: string; x: number; y: number; cluster: string }[];
}

@Component({
  selector: 'app-spatial-correlation-panel',
  standalone: true,
  imports: [
    CommonModule,
    FormsModule,
    MatButtonModule,
    MatButtonToggleModule,
    MatIconModule,
    MatTooltipModule,
    MatProgressSpinnerModule,
  ],
  templateUrl: './spatial-correlation-panel.component.html',
  styleUrls: ['./spatial-correlation-panel.component.scss'],
})
export class SpatialCorrelationPanelComponent implements OnInit, OnChanges {
  @Input() isCompare = false;
  @Input() datasetId?: string;
  /** Whether the app's global split-view compare mode is currently on (not to be confused with
   * `isCompare`, which just says whether THIS panel instance is the compare-side one -- both the
   * main and compare instances get the same value here). Compare mode already shows two real maps
   * side by side, each with its own property, so bivariately recoloring either one would throw
   * that comparison away -- there, the bivariate view instead renders as an independent mini-map
   * inside the Live Correlation drawer. Outside compare mode there's only one map and no drawer,
   * so it's recolored directly instead, like this always did originally. */
  @Input() compareModeActive = false;
  /** Reference geometry (this panel's own dataset's hexagons) for the bivariate mini-map -- this
   * panel doesn't have its own map, it colors the same hexagon shapes the real map uses, just with
   * a 2D blend of featureA/featureB instead of one property. Not wired to cell-click/tab-jump
   * interactions (unlike the real maps) -- kept independent/static for now, same as the drawer's.
   * Named hexFeatures (not `features`) because that name is already the correlatable-properties
   * list below. Only relevant while compareModeActive (see that doc) -- otherwise the real map is
   * recolored directly instead. */
  @Input() hexFeatures: CellFeature[] = [];

  @Output() syncCompareRequested = new EventEmitter<{ featA: SpatialFeature; featB: SpatialFeature }>();
  @Output() applyBivariateColorRequested = new EventEmitter<{ featA: SpatialFeature; featB: SpatialFeature; enabled?: boolean }>();
  @Output() cellSelected = new EventEmitter<string>();

  /** Local UI state for the bivariate toggle button shown outside compare mode (see
   * compareModeActive doc) -- not an @Input(), this button's own on/off state doesn't need to be
   * controlled by the parent the way e.g. datasetId does. */
  isBivariateActive = false;

  @ViewChild('heatmapContainer', { static: false }) heatmapContainer?: ElementRef<HTMLDivElement>;
  @ViewChild('scatterplotContainer', { static: false }) scatterplotContainer?: ElementRef<HTMLDivElement>;
  @ViewChild('bivariateMapContainer', { static: false }) bivariateMapContainer?: ElementRef<HTMLDivElement>;
  @ViewChild('bivariateLegendCanvas', { static: false }) bivariateLegendCanvas?: ElementRef<HTMLCanvasElement>;

  metricMode: 'pearson' | 'moran' = 'pearson';
  isLoadingMatrix = false;
  matrixData: any = null;
  features: SpatialFeature[] = [];

  selectedPairIndices: { i: number; j: number } | null = null;
  selectedFeatureA: SpatialFeature | null = null;
  selectedFeatureB: SpatialFeature | null = null;

  isLoadingScatter = false;
  pairScatterData: PairScatterData | null = null;
  scatterError = '';

  constructor(private http: HttpClient) {}

  ngOnInit(): void {
    if (this.datasetId) {
      this.loadCorrelationMatrix();
    }
  }

  ngOnChanges(changes: SimpleChanges): void {
    if (changes['datasetId'] && this.datasetId) {
      this.loadCorrelationMatrix();
    }
  }

  loadCorrelationMatrix(): void {
    if (!this.datasetId) return;
    this.isLoadingMatrix = true;
    this.matrixData = null;
    this.selectedPairIndices = null;
    this.pairScatterData = null;

    this.http.get<any>(`/api/datasets/${encodeURIComponent(this.datasetId)}/spatial_correlation_matrix`, { withCredentials: true })
      .subscribe({
        next: (res) => {
          this.matrixData = res;
          this.features = res.features || [];
          this.isLoadingMatrix = false;
          setTimeout(() => this.renderHeatmap(), 50);
        },
        error: (err) => {
          console.error('Failed to load spatial correlation matrix:', err);
          this.isLoadingMatrix = false;
        }
      });
  }

  onMetricChange(mode: 'pearson' | 'moran'): void {
    this.metricMode = mode;
    this.renderHeatmap();
  }

  renderHeatmap(): void {
    if (!this.heatmapContainer || !this.matrixData || !this.features.length) return;
    const container = this.heatmapContainer.nativeElement;

    const n = this.features.length;
    if (n === 0) {
      Plotly.purge(container);
      return;
    }

    const xNames = this.features.map(f => f.name);
    const yNames = this.features.map(f => f.name);
    const matrix = this.metricMode === 'pearson' ? this.matrixData.pearson_matrix : this.matrixData.moran_matrix;

    const hoverText: string[][] = [];
    for (let i = 0; i < n; i++) {
      const row: string[] = [];
      for (let j = 0; j < n; j++) {
        const val = matrix[i][j];
        const circ = this.matrixData.circularity_matrix?.[i]?.[j];
        const metricName = this.metricMode === 'pearson' ? 'Pearson r' : "Bivariate Moran's I";
        row.push(
          `<b>${yNames[i]}</b> vs <b>${xNames[j]}</b><br>` +
          `${metricName}: <b>${typeof val === 'number' ? val.toFixed(3) : val}</b><br>` +
          `Status: <i>${circ?.reason || 'Independent'}</i>`
        );
      }
      hoverText.push(row);
    }

    const data: Partial<Plotly.PlotData>[] = [
      {
        type: 'heatmap',
        x: xNames,
        y: yNames,
        z: matrix,
        text: hoverText as any,
        hoverinfo: 'text',
        colorscale: 'RdBu',
        reversescale: true,
        zmin: -1,
        zmax: 1,
        colorbar: {
          title: { text: this.metricMode === 'pearson' ? 'Pearson r' : "Moran's I" },
          thickness: 14,
          len: 0.8
        }
      }
    ];

    const maxNameLen = Math.max(...this.features.map(f => f.name.length), 10);
    const lMargin = Math.max(120, Math.min(220, maxNameLen * 7));
    const bMargin = Math.max(100, Math.min(200, maxNameLen * 6));
    const cellSize = Math.min(32, Math.max(22, Math.floor(500 / n)));
    const matrixSide = n * cellSize;

    const layout: Partial<Plotly.Layout> = {
      width: matrixSide + lMargin + 80,
      height: matrixSide + bMargin + 50,
      margin: {
        t: 40,
        l: lMargin,
        r: 30,
        b: bMargin
      },
      xaxis: {
        tickangle: -45,
        automargin: true,
        type: 'category',
        constrain: 'domain'
      },
      yaxis: {
        automargin: true,
        autorange: 'reversed',
        type: 'category',
        scaleanchor: 'x',
        scaleratio: 1,
        constrain: 'domain'
      }
    };

    Plotly.purge(container);
    Plotly.newPlot(container, data, layout, {
      responsive: true,
      displayModeBar: false
    }).then(() => {
      (container as any).on('plotly_click', (eventData: any) => {
        if (eventData && eventData.points && eventData.points.length > 0) {
          const pt = eventData.points[0];
          const colIdx = pt.x !== undefined ? this.features.findIndex(f => f.name === pt.x) : pt.pointIndex?.[1];
          const rowIdx = pt.y !== undefined ? this.features.findIndex(f => f.name === pt.y) : pt.pointIndex?.[0];
          if (rowIdx >= 0 && colIdx >= 0) {
            this.selectCell(rowIdx, colIdx);
          }
        }
      });
    });
  }

  selectCell(i: number, j: number): void {
    this.selectedPairIndices = { i, j };
    this.selectedFeatureA = this.features[i];
    this.selectedFeatureB = this.features[j];
    this.loadScatterData(this.features[i].id, this.features[j].id);
  }

  loadScatterData(featIdA: string, featIdB: string): void {
    if (!this.datasetId) return;
    this.isLoadingScatter = true;
    this.pairScatterData = null;
    this.scatterError = '';

    const url = `/api/datasets/${encodeURIComponent(this.datasetId)}/spatial_correlation_pair?feature_id_a=${encodeURIComponent(featIdA)}&feature_id_b=${encodeURIComponent(featIdB)}`;

    this.http.get<PairScatterData>(url, { withCredentials: true }).subscribe({
      next: (res) => {
        this.pairScatterData = res;
        this.isLoadingScatter = false;
        setTimeout(() => this.renderScatterplot(), 50);
        // Outside compare mode there's no mini-map to draw -- the real map gets recolored directly
        // instead, via the toggle button (see compareModeActive doc / toggleBivariateMap).
        if (this.compareModeActive) {
          setTimeout(() => this.renderBivariateVisuals(), 50);
        }
      },
      error: (err) => {
        console.error('Failed to load scatter data:', err);
        // Deliberately generic and always the same, regardless of the backend's actual detail --
        // that can be a raw, unfiltered exception message, not something meant for a user to read.
        this.scatterError = "Can't correlate the selected pair";
        this.isLoadingScatter = false;
      }
    });
  }

  renderScatterplot(): void {
    if (!this.scatterplotContainer || !this.pairScatterData) return;
    const scatterData = this.pairScatterData;
    const container = this.scatterplotContainer.nativeElement;

    const data = scatterData.points;
    if (!data || data.length === 0) {
      Plotly.purge(container);
      return;
    }

    const xVals = data.map(d => d.x);
    const yVals = data.map(d => d.y);
    const cellIds = data.map(d => d.cell_id);

    const hoverText = data.map(d =>
      `Cell: <b>${d.cell_id}</b><br>` +
      `Cluster: <b>${d.cluster}</b><br>` +
      `${scatterData.feature_a.name}: <b>${d.x.toFixed(3)}</b><br>` +
      `${scatterData.feature_b.name}: <b>${d.y.toFixed(3)}</b><br>` +
      `<i>Click to highlight on map</i>`
    );

    const colorPalette = d3.schemeCategory10;
    const markerColors = data.map(d => {
      const cNum = parseInt(d.cluster, 10);
      return !isNaN(cNum) ? colorPalette[cNum % 10] : '#1976d2';
    });

    const traces: Partial<Plotly.PlotData>[] = [
      {
        type: 'scatter',
        mode: 'markers',
        x: xVals,
        y: yVals,
        text: hoverText,
        hoverinfo: 'text',
        customdata: cellIds as any,
        marker: {
          size: 6,
          color: markerColors,
          opacity: 0.65
        },
        name: 'Cells'
      }
    ];

    // Linear Regression Trendline
    const { slope, intercept } = scatterData.stats;
    if (slope !== 0 || intercept !== 0) {
      const minX = Math.min(...xVals);
      const maxX = Math.max(...xVals);
      traces.push({
        type: 'scatter',
        mode: 'lines',
        x: [minX, maxX],
        y: [slope * minX + intercept, slope * maxX + intercept],
        line: { color: '#e53935', width: 2, dash: 'dash' },
        name: 'OLS Trend',
        hoverinfo: 'none'
      });
    }

    const layout: Partial<Plotly.Layout> = {
      margin: { t: 20, l: 60, r: 20, b: 60 },
      height: 320,
      xaxis: {
        title: { text: scatterData.feature_a.name },
        automargin: true
      },
      yaxis: {
        title: { text: scatterData.feature_b.name },
        automargin: true
      },
      showlegend: false,
      hovermode: 'closest'
    };

    Plotly.purge(container);
    Plotly.newPlot(container, traces, layout, {
      responsive: true,
      displayModeBar: false
    }).then(() => {
      (container as any).on('plotly_click', (eventData: any) => {
        if (eventData && eventData.points && eventData.points.length > 0) {
          const pt = eventData.points[0];
          const cellId = pt.customdata;
          if (cellId) {
            this.cellSelected.emit(cellId);
          }
        }
      });
    });
  }

  onSyncCompare(): void {
    if (this.selectedFeatureA && this.selectedFeatureB) {
      this.syncCompareRequested.emit({
        featA: this.selectedFeatureA,
        featB: this.selectedFeatureB
      });
    }
  }

  /** Outside compare mode only (see compareModeActive doc) -- recolors the real map directly via
   * the parent's onApplyBivariateColor. */
  toggleBivariateMap(): void {
    if (this.selectedFeatureA && this.selectedFeatureB) {
      this.isBivariateActive = !this.isBivariateActive;
      this.applyBivariateColorRequested.emit({
        featA: this.selectedFeatureA,
        featB: this.selectedFeatureB,
        enabled: this.isBivariateActive
      });
    }
  }

  /** Same formula as live-correlation-drawer.component.ts's bivariateColor -- kept in sync by
   * hand since these are two separate components with their own render pipelines (matching this
   * codebase's existing pattern of each having its own renderScatterplot rather than sharing one).
   * u/v are featureA/featureB each normalized to [0,1]; green marks how close they are to each
   * other. Corners: (0,0)=green, (1,1)=white, (0,1)=red, (1,0)=blue. */
  private bivariateColor(u: number, v: number): string {
    const r = Math.round(255 * v);
    const g = Math.round(255 * Math.max(0, 1 - Math.abs(u - v)));
    const b = Math.round(255 * u);
    return `rgb(${r},${g},${b})`;
  }

  private renderBivariateVisuals(): void {
    this.renderBivariateMap();
    this.renderBivariateLegend();
  }

  /** Renders hexFeatures (this panel's own dataset's real hexagon shapes) colored by a 2D blend of
   * featureA (x) / featureB (y) instead of a single property. Uses bivariate_coords, which the
   * same /spatial_correlation_pair fetch behind the scatterplot above already returns for every
   * cell (not just the scatter sample), so this needs no separate request. Static/non-interactive
   * for now -- no click-to-select, no tab-jump wiring, matching the real maps. */
  private renderBivariateMap(): void {
    const container = this.bivariateMapContainer?.nativeElement;
    const coords = (this.pairScatterData as any)?.bivariate_coords;
    if (!container || !coords?.cell_ids?.length || !this.hexFeatures?.length) return;

    const pointsMap = new Map<string, { x: number; y: number }>();
    for (let i = 0; i < coords.cell_ids.length; i++) {
      pointsMap.set(String(coords.cell_ids[i]), { x: coords.x[i], y: coords.y[i] });
    }
    const allX: number[] = coords.x;
    const allY: number[] = coords.y;
    const minX = Math.min(...allX), maxX = Math.max(...allX);
    const minY = Math.min(...allY), maxY = Math.max(...allY);

    const width = container.clientWidth || 300;
    const height = width;
    const padding = 10;

    d3.select(container).selectAll('svg').remove();
    const svg = d3.select(container)
      .append('svg')
      .attr('width', width)
      .attr('height', height)
      .attr('viewBox', `0 0 ${width} ${height}`);

    const projection = d3.geoIdentity().fitExtent(
      [[padding, padding], [width - padding, height - padding]],
      { type: 'FeatureCollection', features: this.hexFeatures } as any,
    );
    const pathGenerator = d3.geoPath<CellFeature>().projection(projection);

    svg.selectAll('path')
      .data(this.hexFeatures, (d: any) => d.properties.barcode)
      .join('path')
      .attr('d', (d) => pathGenerator(d) || '')
      .attr('fill', (d) => {
        const cellId = String(d.properties?.barcode ?? d.id ?? '');
        const pt = pointsMap.get(cellId);
        if (!pt) return '#ccc';
        const u = maxX > minX ? (pt.x - minX) / (maxX - minX) : 0.5;
        const v = maxY > minY ? (pt.y - minY) / (maxY - minY) : 0.5;
        return this.bivariateColor(u, v);
      })
      .style('stroke', 'none')
      .style('opacity', 0.9);
  }

  /** Draws the actual 2D blend as a pixel grid (bivariateColor isn't a simple 2-stop gradient, so
   * no CSS linear-gradient can represent it) so this legend always matches renderBivariateMap
   * exactly instead of being a separately hand-picked gradient that can drift out of sync. */
  private renderBivariateLegend(): void {
    const canvas = this.bivariateLegendCanvas?.nativeElement;
    if (!canvas) return;
    const size = 64;
    canvas.width = size;
    canvas.height = size;
    const ctx = canvas.getContext('2d');
    if (!ctx) return;

    const image = ctx.createImageData(size, size);
    for (let row = 0; row < size; row++) {
      const v = 1 - row / (size - 1);
      for (let col = 0; col < size; col++) {
        const u = col / (size - 1);
        const [r, g, b] = this.bivariateColor(u, v)
          .replace(/^rgb\(|\)$/g, '')
          .split(',')
          .map(Number);
        const idx = (row * size + col) * 4;
        image.data[idx] = r;
        image.data[idx + 1] = g;
        image.data[idx + 2] = b;
        image.data[idx + 3] = 255;
      }
    }
    ctx.putImageData(image, 0, 0);
  }
}
