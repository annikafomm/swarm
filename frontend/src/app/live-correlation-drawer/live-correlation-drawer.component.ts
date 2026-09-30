import {
  Component,
  Input,
  Output,
  EventEmitter,
  OnChanges,
  SimpleChanges,
  ViewChild,
  ElementRef,
  OnDestroy
} from '@angular/core';
import { CommonModule } from '@angular/common';
import { HttpClient } from '@angular/common/http';
import { MatIconModule } from '@angular/material/icon';
import { MatButtonModule } from '@angular/material/button';
import { MatTooltipModule } from '@angular/material/tooltip';
import { MatProgressSpinnerModule } from '@angular/material/progress-spinner';
import * as Plotly from 'plotly.js-dist-min';
import * as d3 from 'd3';
import { CellFeature } from '../hexagon-view/cell-feature.types';

export interface SpatialFeatureRef {
  id: string;
  name: string;
}

@Component({
  selector: 'app-live-correlation-drawer',
  standalone: true,
  imports: [
    CommonModule,
    MatIconModule,
    MatButtonModule,
    MatTooltipModule,
    MatProgressSpinnerModule
  ],
  templateUrl: './live-correlation-drawer.component.html',
  styleUrls: ['./live-correlation-drawer.component.scss']
})
export class LiveCorrelationDrawerComponent implements OnChanges, OnDestroy {
  @Input() datasetId: string | undefined;
  /** The compare map's dataset, when featureB comes from a different dataset than featureA
   * (e.g. main vs. compare showing separate datasets rather than two properties of one).
   * Defaults to `datasetId` when not set, matching the same-dataset case. */
  @Input() datasetIdCompare: string | undefined;
  @Input() featureA: SpatialFeatureRef | null = null;
  @Input() featureB: SpatialFeatureRef | null = null;
  @Input() isVisible: boolean = false;
  /** Reference geometry (main dataset's hexagons) for the bivariate mini-map -- this drawer
   * doesn't have its own map, it colors the SAME hexagon shapes the real main map uses, just with
   * a 2D blend of featureA/featureB instead of one property. Not wired to cell-click/tab-jump
   * interactions (unlike the real maps) -- kept independent/static for now. */
  @Input() features: CellFeature[] = [];

  @Output() cellSelected = new EventEmitter<string>();
  @Output() drawerClosed = new EventEmitter<void>();

  @ViewChild('drawerScatterContainer') scatterContainer!: ElementRef<HTMLDivElement>;
  @ViewChild('bivariateMapContainer') bivariateMapContainer?: ElementRef<HTMLDivElement>;
  @ViewChild('bivariateLegendCanvas') bivariateLegendCanvas?: ElementRef<HTMLCanvasElement>;

  public isCollapsed: boolean = false;
  public isLoading: boolean = false;
  public scatterData: any = null;
  public errorMsg: string = '';

  private currentFetchSub: any = null;

  constructor(private http: HttpClient) {}

  public onClose(): void {
    this.drawerClosed.emit();
  }

  ngOnChanges(changes: SimpleChanges): void {
    const prevA = changes['featureA']?.previousValue?.id;
    const currA = changes['featureA']?.currentValue?.id;
    const prevB = changes['featureB']?.previousValue?.id;
    const currB = changes['featureB']?.currentValue?.id;
    const prevDs = changes['datasetId']?.previousValue;
    const currDs = changes['datasetId']?.currentValue;
    const prevDsB = changes['datasetIdCompare']?.previousValue;
    const currDsB = changes['datasetIdCompare']?.currentValue;

    const featAChanged = prevA !== currA && !!currA;
    const featBChanged = prevB !== currB && !!currB;
    const datasetChanged = (prevDs !== currDs && !!currDs) || prevDsB !== currDsB;
    const visChanged = changes['isVisible'] && changes['isVisible'].currentValue && !changes['isVisible'].previousValue;

    if (this.isVisible && (featAChanged || featBChanged || datasetChanged || visChanged)) {
      this.loadCorrelationData();
    }
  }

  ngOnDestroy(): void {
    if (this.scatterContainer?.nativeElement) {
      Plotly.purge(this.scatterContainer.nativeElement);
    }
  }

  public toggleCollapse(): void {
    this.isCollapsed = !this.isCollapsed;
    if (!this.isCollapsed && this.scatterData) {
      setTimeout(() => this.renderScatterplot(), 50);
      setTimeout(() => this.renderBivariateVisuals(), 50);
    }
  }

  public loadCorrelationData(): void {
    if (!this.datasetId || !this.featureA?.id || !this.featureB?.id) {
      return;
    }

    this.isLoading = true;
    this.errorMsg = '';

    let url = `/api/datasets/${encodeURIComponent(this.datasetId)}/spatial_correlation_pair?feature_id_a=${encodeURIComponent(this.featureA.id)}&feature_id_b=${encodeURIComponent(this.featureB.id)}`;
    if (this.datasetIdCompare && this.datasetIdCompare !== this.datasetId) {
      url += `&dataset_id_b=${encodeURIComponent(this.datasetIdCompare)}`;
    }

    if (this.currentFetchSub) {
      this.currentFetchSub.unsubscribe();
    }

    this.currentFetchSub = this.http.get<any>(url, { withCredentials: true }).subscribe({
      next: (res) => {
        this.isLoading = false;
        this.scatterData = res;
        if (!this.isCollapsed) {
          setTimeout(() => this.renderScatterplot(), 50);
        }
        setTimeout(() => this.renderBivariateVisuals(), 50);
      },
      error: (err) => {
        this.isLoading = false;
        // Deliberately generic and always the same message, regardless of the backend's actual
        // detail -- that detail can be a raw, unfiltered exception message (e.g. an internal
        // TypeError from a computation bug), not something meant for an end user to read.
        this.errorMsg = "Can't correlate the selected pair";
        console.error('Failed to load drawer scatter data:', err);
      }
    });
  }

  private renderScatterplot(): void {
    if (!this.scatterContainer?.nativeElement || !this.scatterData) return;

    const container = this.scatterContainer.nativeElement;
    const points = this.scatterData.points || [];
    const stats = this.scatterData.stats || {};

    const xVals = points.map((p: any) => p.x);
    const yVals = points.map((p: any) => p.y);
    const cellIds = points.map((p: any) => p.cell_id);
    const clusters = points.map((p: any) => p.cluster);

    // `cluster` is a label (leiden id or cell_type name, both stringified by the backend), not a
    // number -- passing it straight into `color` with a continuous `colorscale` (as before) only
    // works by accident for numeric-looking leiden ids, and silently falls back to a single flat
    // color for cell_type labels like "B cell", since those aren't valid CSS color names either.
    // Map each distinct cluster to its own color from a qualitative palette instead.
    const clusterColorScale = d3.scaleOrdinal<string>(d3.schemeSet2).domain(Array.from(new Set(clusters)));
    const pointColors = clusters.map((c: string) => clusterColorScale(c));

    const hoverTexts = points.map(
      (p: any) =>
        `<b>Cell:</b> ${p.cell_id}<br>` +
        `<b>${this.featureA?.name}:</b> ${typeof p.x === 'number' ? p.x.toFixed(3) : p.x}<br>` +
        `<b>${this.featureB?.name}:</b> ${typeof p.y === 'number' ? p.y.toFixed(3) : p.y}<br>` +
        `<b>Cluster:</b> ${p.cluster}`
    );

    const data: Partial<Plotly.PlotData>[] = [
      {
        type: 'scatter',
        mode: 'markers',
        x: xVals,
        y: yVals,
        text: hoverTexts,
        customdata: cellIds,
        hoverinfo: 'text',
        marker: {
          size: 5,
          color: pointColors,
          opacity: 0.75
        },
        name: 'Cells'
      }
    ];

    if (points.length > 1 && stats.slope !== undefined && stats.intercept !== undefined) {
      const minX = Math.min(...xVals);
      const maxX = Math.max(...xVals);
      data.push({
        type: 'scatter',
        mode: 'lines',
        x: [minX, maxX],
        y: [minX * stats.slope + stats.intercept, maxX * stats.slope + stats.intercept],
        line: { color: '#ef4444', width: 2, dash: 'dot' },
        name: 'OLS Trend'
      });
    }

    const layout: Partial<Plotly.Layout> = {
      margin: { t: 15, l: 45, r: 15, b: 35 },
      height: 220,
      showlegend: false,
      xaxis: {
        title: { text: this.featureA?.name || 'Map 1 (Left)', font: { size: 10 } },
        automargin: true,
        tickfont: { size: 9 }
      },
      yaxis: {
        title: { text: this.featureB?.name || 'Map 2 (Right)', font: { size: 10 } },
        automargin: true,
        tickfont: { size: 9 }
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
          if (pt.customdata) {
            this.cellSelected.emit(String(pt.customdata));
          }
        }
      });
    });
  }

  /** The one place the featureA(u)/featureB(v) -> color blend is defined -- both the mini-map and
   * the legend canvas call this, so the legend is guaranteed to actually match what the map shows
   * instead of being a separately hand-picked gradient that can drift out of sync with it (the
   * previous static CSS gradient swatch was a blue->magenta look that didn't match this formula's
   * actual corners at all: (0,0)=green, (1,1)=white, (0,1)=red, (1,0)=blue). u and v are each
   * featureA/featureB normalized to [0,1]; green marks how close the two are to each other. */
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

  /** Renders the hexagon geometry from `features` (the main map's own shapes -- this drawer has
   * no map of its own) colored by a 2D blend of featureA (x) / featureB (y) instead of a single
   * property. Uses `bivariate_coords`, which the same /spatial_correlation_pair fetch behind the
   * scatterplot already returns for every cell (not just the <=2000-point scatter sample), so this
   * needs no separate network request. Static/non-interactive for now -- no click-to-select, no
   * tab-jump wiring, matching the real maps' interactions. */
  private renderBivariateMap(): void {
    const container = this.bivariateMapContainer?.nativeElement;
    const coords = this.scatterData?.bivariate_coords;
    if (!container || !coords?.cell_ids?.length || !this.features?.length) return;

    const pointsMap = new Map<string, { x: number; y: number }>();
    for (let i = 0; i < coords.cell_ids.length; i++) {
      pointsMap.set(String(coords.cell_ids[i]), { x: coords.x[i], y: coords.y[i] });
    }
    const allX: number[] = coords.x;
    const allY: number[] = coords.y;
    const minX = Math.min(...allX), maxX = Math.max(...allX);
    const minY = Math.min(...allY), maxY = Math.max(...allY);

    const width = container.clientWidth || 300;
    const height = width; // square, matching the real map's own aspect ratio convention
    const padding = 10;

    d3.select(container).selectAll('svg').remove();
    const svg = d3.select(container)
      .append('svg')
      .attr('width', width)
      .attr('height', height)
      .attr('viewBox', `0 0 ${width} ${height}`);

    const projection = d3.geoIdentity().fitExtent(
      [[padding, padding], [width - padding, height - padding]],
      { type: 'FeatureCollection', features: this.features } as any,
    );
    const pathGenerator = d3.geoPath<CellFeature>().projection(projection);

    svg.selectAll('path')
      .data(this.features, (d: any) => d.properties.barcode)
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

  /** Draws the actual 2D blend as a pixel grid (bivariateColor isn't a simple 2-stop gradient --
   * the green channel depends on |u-v|, so no CSS linear-gradient can represent it correctly) so
   * the legend always matches renderBivariateMap exactly, corner for corner. */
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
      // Canvas y grows downward; v (featureB) should read bottom-to-top like the axis label does.
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

