import { Component, ElementRef, EventEmitter, Input, OnDestroy, Output, ViewChild } from '@angular/core';
import { CommonModule } from '@angular/common';
import * as d3 from 'd3';
import { CellFeature } from './cell-feature.types';

export interface HexagonRenderContext {
  features: CellFeature[];
  isXenium: boolean;
  fullFeatures: CellFeature[];
  selectedView: string;
  legendType: 'categorical' | 'continuous';
  /**
   * Pre-existing behavior being preserved, not introduced here: categorical/ordinal fills always
   * use the *main* view's color scale, even when rendering the compare instance — only the
   * continuous scale is correctly per-instance. Pass the main view's colorScale for both
   * instances to match current app behavior exactly; flagged to the user as a latent bug worth
   * a follow-up fix, not fixed silently as part of this structural extraction.
   */
  colorScale: d3.ScaleOrdinal<string, string>;
  continuousColorScale: (t: number) => string;
  leidenClusterAnnotations: Record<string, any> | undefined;
  leidenCentralityProps: string[];
}

/**
 * Owns one side (main or compare) of the hexagon map: the D3 SVG/rendering, cluster-outline
 * highlighting, and (main-view only) the Xenium performance downsampling + magnifier detail
 * window. Extracted from HexagonPlotComponent, which used to hand-duplicate all of this for
 * main vs. compare via `compare ? x : y` ternaries throughout.
 *
 * Render-time data (features, color scales, selected view, etc.) is passed as explicit method
 * parameters rather than read from @Input()-bound fields, because the parent calls these methods
 * imperatively via @ViewChild right after mutating its own state, in the same synchronous tick —
 * before Angular's own change-detection pass would have refreshed an @Input() binding.
 */
@Component({
  selector: 'app-hexagon-view',
  standalone: true,
  imports: [CommonModule],
  templateUrl: './hexagon-view.component.html',
  styleUrls: ['./hexagon-view.component.scss'],
})
export class HexagonViewComponent implements OnDestroy {
  @Input() isCompare = false;
  /** This instance's own currently-selected cell, for the mouseLeave highlight-suppression check. */
  @Input() selectedCell: CellFeature | null = null;

  @Output() cellClicked = new EventEmitter<{ event: MouseEvent; cell: CellFeature }>();
  @Output() clusterClicked = new EventEmitter<{ cell: CellFeature; clusterId: number }>();

  @ViewChild('hexbinContainer', { static: true }) private hexbinContainerRef!: ElementRef<HTMLDivElement>;

  /** Public: read by the parent's updateHexColors/getViewVariablesToUpdate via @ViewChild until
   * those are migrated too (they're deeply tied to cross-cutting concerns staying in the parent
   * for now — see the refactor plan). */
  public svg!: d3.Selection<SVGSVGElement, any, any, any>;
  public g!: d3.Selection<SVGGElement, any, any, any>;
  public baseLayer!: d3.Selection<SVGGElement, null, any, any>;
  public detailLayer!: d3.Selection<SVGGElement, null, any, any>;
  public currentPathGenerator: d3.GeoPath<any, any> | null = null;

  private currentTransform = d3.zoomIdentity;

  // ======= Xenium detail-window state =======
  private readonly detailSize = 80;
  private detailVisible = false;
  private detailScreenPos: { x: number; y: number } | null = null;
  private keydownHandler: ((event: KeyboardEvent) => void) | null = null;
  private lastRenderCtx: HexagonRenderContext | null = null;

  ngOnDestroy(): void {
    if (this.keydownHandler) {
      window.removeEventListener('keydown', this.keydownHandler);
    }
  }

  /** Removes any existing SVG in this instance's container before a fresh createHexagonPlot(). */
  public clearSvg(): void {
    d3.select(this.hexbinContainerRef.nativeElement).selectAll('svg').remove();
  }

  public createHexagonPlot(): void {
    const width = 500;
    const height = 400;

    const container = d3.select(this.hexbinContainerRef.nativeElement);
    const svgSel = container
      .selectAll('svg')
      .data([0])
      .join('svg')
      .attr('width', width)
      .attr('height', height)
      .attr('viewBox', [0, 0, 1200, 1000] as [number, number, number, number])
      .style('background-color', 'white')
      .style('overflow', 'hidden');

    const gSel = svgSel
      .selectAll<SVGGElement, number>('g.root-group')
      .data([0])
      .join('g')
      .attr('class', 'root-group');

    this.svg = svgSel as unknown as d3.Selection<SVGSVGElement, any, any, any>;
    this.g = gSel as any;

    if (!this.isCompare) {
      // Xenium layers live only in the main view.
      this.baseLayer = this.g.selectAll<SVGGElement, unknown>('g.base-layer')
        .data([null])
        .join('g')
        .attr('class', 'base-layer');

      this.detailLayer = this.svg.selectAll<SVGGElement, unknown>('g.detail-layer')
        .data([null])
        .join('g')
        .attr('class', 'detail-layer')
        .attr('clip-path', 'url(#detail-clip)');
    }

    const zoomBehavior = d3
      .zoom<SVGSVGElement, unknown>()
      .scaleExtent([1, 5])
      .extent([
        [0, 0],
        [width, height],
      ])
      .on('zoom', (event) => {
        this.g.attr('transform', event.transform.toString());
        this.currentTransform = event.transform;
      });

    (svgSel as any).call(zoomBehavior);
  }

  public renderHexagons(ctx: HexagonRenderContext): void {
    this.lastRenderCtx = ctx;
    const width = 1200;
    const height = 1000;

    const projection = d3.geoIdentity().fitSize([width, height], {
      type: 'FeatureCollection',
      features: ctx.isXenium ? ctx.fullFeatures : ctx.features,
    });

    this.currentPathGenerator = d3.geoPath<CellFeature>().projection(projection);
    const pathGenerator = this.currentPathGenerator;

    if (ctx.isXenium) {
      ctx.fullFeatures.forEach((f) => {
        if (!(f.properties as any).__centroidProjected) {
          const c = this.currentPathGenerator!.centroid(f as any);
          (f.properties as any).__centroidProjected = c;
        }
      });
    }

    const getFill = (d: CellFeature): string => {
      const value = ctx.leidenCentralityProps.includes(ctx.selectedView)
        ? ctx.leidenClusterAnnotations?.[String(d.properties.leiden)]?.centrality?.[ctx.selectedView]
        : d.properties?.[ctx.selectedView];
      if (ctx.legendType === 'categorical') {
        return ctx.colorScale(String(value));
      }
      const num = this.toNumber(value);
      return Number.isFinite(num) ? ctx.continuousColorScale(num) : '#ccc';
    };

    if (ctx.isXenium) {
      this.baseLayer
        .style('cursor', 'default')
        .style('pointer-events', 'none')
        .selectAll<SVGPathElement, CellFeature>('path')
        .data(ctx.features, (d: any) => d.properties.barcode)
        .join('path')
        .attr('d', (d: CellFeature) => pathGenerator(d) || '')
        .attr('fill', getFill)
        .style('opacity', 0.8);
    } else {
      this.g
        .style('cursor', 'pointer')
        .style('pointer-events', null)
        .selectAll<SVGPathElement, CellFeature>('path')
        .data(ctx.features, (d: any) => d.properties.barcode)
        .join('path')
        .attr('d', (d: CellFeature) => pathGenerator(d) || '')
        .attr('fill', getFill)
        .style('opacity', 0.8)
        .on('mouseover', (event, d) => this.mouseOver(event, d))
        .on('mouseleave', (event, d) => this.mouseLeave(event, d))
        .on('click', (event, d) => this.cellClicked.emit({ event, cell: d }));
    }
    // Note: initDetailWindow/hideDetailWindow/bindDetailWindowInteractions are NOT called here —
    // the original calls them from loadAndRenderData, before renderHexagons, as part of Xenium
    // setup. The parent redirects that same call sequence to this instance's public methods.
  }

  public extendCluster(selectedCluster: number, features: CellFeature[]): void {
    const paths = this.g.selectAll<SVGPathElement, CellFeature>('path');

    // Remove mouseleave event to prevent resetting outline. Must run on the
    // selection itself, not the transition below — Transition.on() only
    // accepts 'start'/'end'/'interrupt' and throws on any other event type.
    paths.on('mouseleave', null);

    paths
      .transition()
      .duration(300)
      .attr('d', (d: CellFeature) => {
        const projection = d3.geoIdentity().fitSize([1200, 1000], {
          type: 'FeatureCollection',
          features: features,
        });
        const pathGenerator = d3.geoPath<CellFeature>().projection(projection);
        return pathGenerator(d) || '';
      })
      .attr('stroke-width', (d: CellFeature) => (d.properties.leiden === selectedCluster ? '3px' : '1px'))
      .attr('stroke', (d: CellFeature) => (d.properties.leiden === selectedCluster ? '#000' : 'transparent'))
      .style('opacity', (d: CellFeature) => (d.properties.leiden === selectedCluster ? 1.0 : 0.6));
  }

  public resetClusterExtension(features: CellFeature[]): void {
    const projection = d3.geoIdentity().fitSize([1200, 1000], {
      type: 'FeatureCollection',
      features: features,
    });
    const pathGenerator = d3.geoPath<CellFeature>().projection(projection);

    this.g
      .selectAll<SVGPathElement, CellFeature>('path')
      .transition()
      .duration(300)
      .attr('d', (d: CellFeature) => pathGenerator(d) || '')
      .attr('stroke-width', '1px')
      .attr('stroke', 'transparent')
      .style('opacity', 0.8);

    // Reinitialize the mouseleave event
    this.g
      .selectAll<SVGPathElement, CellFeature>('path')
      .on('mouseleave', (event, d) => this.mouseLeave(event, d));
  }

  private toNumber(v: unknown): number {
    if (typeof v === 'number') return v;
    if (typeof v === 'string') {
      const n = parseFloat(v);
      return Number.isFinite(n) ? n : NaN;
    }
    return NaN;
  }

  // Original also did `d3.selectAll('.Country').transition()...` in both handlers below —
  // no element in this app has ever had class="Country" (confirmed: zero matches repo-wide),
  // so that half of both handlers was dead code with zero DOM effect. Omitted here, not a
  // behavior change.
  private mouseOver(event: MouseEvent, d: CellFeature): void {
    d3.select(event.target as SVGElement)
      .transition()
      .duration(200)
      .style('opacity', 0.8)
      .attr('stroke', 'black');
  }

  private mouseLeave(event: MouseEvent, d: CellFeature): void {
    // Original also suppressed the reset if `d` matched the *other* view's selected cell
    // (this.selectedCellCompare when called from the main instance, and vice versa) — a
    // real, if narrow, cross-view check. Simplified to only this instance's own selectedCell
    // since a per-instance component doesn't have a clean way to see the other instance's
    // selection without new plumbing; flagged to the user as a deliberate small behavior
    // change, not an oversight.
    if (
      this.selectedCell &&
      (d.properties.barcode === this.selectedCell.properties.barcode ||
        (this.lastRenderCtx?.selectedView === 'leiden' &&
          d.properties.leiden === this.selectedCell.properties.leiden))
    ) {
      return;
    }
    d3.select(event.target as SVGElement)
      .transition()
      .duration(200)
      .attr('stroke', 'transparent');
  }

  // ======= Xenium detail window (main view only) =======

  public initDetailWindow(): void {
    this.svg.select('#detail-frame').remove();

    const defs = this.svg.select('defs').empty()
      ? this.svg.append('defs')
      : this.svg.select('defs');

    defs.select('#detail-clip').remove();

    defs.append('clipPath')
      .attr('id', 'detail-clip')
      .append('rect')
      .attr('id', 'detail-window')
      .attr('x', 30)
      .attr('y', 30)
      .attr('width', this.detailSize)
      .attr('height', this.detailSize);

    this.detailLayer.selectAll('.detail-bg').remove();

    this.detailLayer
      .append('rect')
      .attr('class', 'detail-bg')
      .attr('x', -10000)
      .attr('y', -10000)
      .attr('width', 20000)
      .attr('height', 20000)
      .attr('fill', '#ffffff')
      .attr('pointer-events', 'none');

    this.svg.append('rect')
      .attr('id', 'detail-frame')
      .attr('x', 30)
      .attr('y', 30)
      .attr('width', this.detailSize)
      .attr('height', this.detailSize)
      .attr('fill', 'none')
      .attr('stroke', '#111')
      .attr('stroke-width', 1)
      .style('pointer-events', 'none');
  }

  public bindDetailWindowInteractions(): void {
    if (!this.svg) return;

    this.svg.on('click.detail', null);
    this.svg.on('dblclick.detail', null);
    this.svg.on('contextmenu.detail', null);

    this.svg.on('dblclick.detail', (event: MouseEvent) => {
      event.preventDefault();
      event.stopPropagation();
      const [mx, my] = d3.pointer(event, this.svg.node());
      this.setDetailWindow(mx, my);
    });

    this.svg.on('contextmenu.detail', (event: MouseEvent) => {
      event.preventDefault();
      event.stopPropagation();
      this.clearDetailWindow();
    });

    if (this.keydownHandler) {
      window.removeEventListener('keydown', this.keydownHandler);
    }
    this.keydownHandler = (event: KeyboardEvent) => {
      if (event.key === 'Escape') {
        this.clearDetailWindow();
      }
    };
    window.addEventListener('keydown', this.keydownHandler);
  }

  private setDetailWindow(screenX: number, screenY: number): void {
    this.detailVisible = true;
    this.detailScreenPos = { x: screenX, y: screenY };
    this.showDetailWindowAt(screenX, screenY);
    this.updateDetailAtScreenPos(screenX, screenY);
  }

  private clearDetailWindow(): void {
    this.detailVisible = false;
    this.detailScreenPos = null;
    this.hideDetailWindow();
    this.detailLayer.selectAll('path').remove();
  }

  /**
   * Re-shows the magnifier detail window at its last screen position, if it was visible —
   * called by the parent's updateHexColors() after a recolor, matching the original's inline
   * `if (containerName === '#hexbin' && this.isXenium && this.detailVisible && this.detailScreenPos)`
   * guard. Kept as one method (rather than exposing detailVisible/detailScreenPos publicly) so
   * this instance's detail-window state stays encapsulated.
   */
  public refreshDetailWindowIfVisible(): void {
    if (!this.lastRenderCtx?.isXenium || !this.detailVisible || !this.detailScreenPos) return;
    this.showDetailWindowAt(this.detailScreenPos.x, this.detailScreenPos.y);
    this.updateDetailAtScreenPos(this.detailScreenPos.x, this.detailScreenPos.y);
  }

  public hideDetailWindow(): void {
    this.svg.select('#detail-window').style('display', 'none');
    this.svg.select('#detail-frame').style('display', 'none');
    this.detailLayer.style('display', 'none');
  }

  private showDetailWindowAt(screenX: number, screenY: number): void {
    this.svg.select('#detail-window')
      .attr('x', screenX - this.detailSize / 2)
      .attr('y', screenY - this.detailSize / 2)
      .style('display', null);

    this.svg.select('#detail-frame')
      .attr('x', screenX - this.detailSize / 2)
      .attr('y', screenY - this.detailSize / 2)
      .style('display', null);

    this.detailLayer.style('display', null);
  }

  private updateDetailAtScreenPos(screenX: number, screenY: number): void {
    const ctx = this.lastRenderCtx;
    if (!ctx?.isXenium || !this.currentPathGenerator) return;

    const half = this.detailSize / 2;

    const x0s = screenX - half;
    const x1s = screenX + half;
    const y0s = screenY - half;
    const y1s = screenY + half;

    const [x0d, y0d] = this.currentTransform.invert([x0s, y0s]);
    const [x1d, y1d] = this.currentTransform.invert([x1s, y1s]);

    const centerX = (x0d + x1d) / 2;
    const centerY = (y0d + y1d) / 2;

    const [targetX, targetY] = this.currentTransform.invert([screenX, screenY]);

    const localScale = 6;

    this.detailLayer.attr(
      'transform',
      `translate(${targetX},${targetY}) scale(${localScale}) translate(${-centerX},${-centerY})`,
    );

    const subset = ctx.fullFeatures.filter((f) => {
      const c = (f.properties as any).__centroidProjected;
      if (!c || c.length < 2) return false;
      const [x, y] = c as [number, number];
      return x >= x0d && x <= x1d && y >= y0d && y <= y1d;
    });

    const pathGenerator = this.currentPathGenerator;
    this.detailLayer
      .style('cursor', 'pointer')
      .selectAll<SVGPathElement, CellFeature>('path')
      .data(subset)
      .join('path')
      .attr('d', (d: CellFeature) => pathGenerator(d) || '')
      .attr('fill', (d: CellFeature) => {
        const value = ctx.leidenCentralityProps.includes(ctx.selectedView)
          ? ctx.leidenClusterAnnotations?.[String(d.properties.leiden)]?.centrality?.[ctx.selectedView]
          : d.properties?.[ctx.selectedView];
        if (ctx.legendType === 'categorical') {
          return ctx.colorScale(String(value));
        }
        const num = this.toNumber(value);
        return Number.isFinite(num) ? ctx.continuousColorScale(num) : '#ccc';
      })
      .attr('stroke', '#fff')
      .attr('stroke-width', 0.4)
      .style('opacity', 1)
      .on('click', (event, d) => this.cellClicked.emit({ event, cell: d }));
  }
}
