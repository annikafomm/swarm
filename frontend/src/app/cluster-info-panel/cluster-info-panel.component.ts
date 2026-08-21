import { Component, EventEmitter, Input, Output } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import * as Plotly from 'plotly.js-dist-min';

export interface ClusterCellTypeCount {
  type: string;
  count: number;
  percentage: string;
}

export interface ClusterCentralityAvg {
  degree_centrality: number;
  average_clustering: number;
  closeness_centrality: number;
}

/**
 * Renders the "Cluster Information" sidebar tab for one side (main or compare) — cluster
 * stats, cell-type distribution, centrality averages, and the neighborhood-enrichment heatmap.
 * Extracted from HexagonPlotComponent, which used to hand-duplicate this markup verbatim for
 * main vs. compare (no ngTemplateOutlet sharing, unlike the Cell Info panel).
 *
 * Unlike HexagonViewComponent/CellInfoPanelComponent, this is a genuinely passive, @Input()-
 * driven display — the actual cluster-selection orchestration (displayClusterDetails,
 * selectCluster, calculateClusterStats, clearClusterData, autoSelectDefaultCluster) stays in
 * HexagonPlotComponent, because it also touches co-occurrence/GRN/footprint state that belongs
 * to later refactor phases. Only renderNhoodHeatmap (self-contained: reads its inputs, draws
 * into its own DOM id) moves here, called by the parent's onTabChange via @ViewChild.
 */
import { MatIconModule } from '@angular/material/icon';
import { MatButtonModule } from '@angular/material/button';
import { MatTooltipModule } from '@angular/material/tooltip';
import { InfoService } from '../info.service';

@Component({
  selector: 'app-cluster-info-panel',
  standalone: true,
  imports: [CommonModule, FormsModule, MatIconModule, MatButtonModule, MatTooltipModule],
  templateUrl: './cluster-info-panel.component.html',
  styleUrls: ['./cluster-info-panel.component.scss'],
})
export class ClusterInfoPanelComponent {
  constructor(public infoService: InfoService) {}

  @Input() isCompare = false;
  /** Mirrors the original's extra `selectedView === 'leiden'` guard on top of "a cluster tab
   * exists at all" — the tab can exist while the map is colored by something else. When false
   * (and a cluster is selected), the template shows an explanatory message instead of the full
   * stats, since none of the stats below apply outside the leiden coloring. */
  @Input() isLeidenView = false;
  @Input() clusterId: number | null = null;
  /** All leiden cluster ids for this dataset, to populate the cluster-select dropdown. */
  @Input() availableClusterIds: number[] = [];
  /**
   * Emits when the user picks a different cluster from the dropdown. The parent is expected to
   * call its own selectCluster(id, isCompare) in response — not just reassign clusterId — since
   * that also updates clusterCells/co-occurrence/map highlight; this component stays passive.
   */
  @Output() clusterSelected = new EventEmitter<number>();
  @Input() clusterCellCount = 0;
  @Input() clusterCellTypes: ClusterCellTypeCount[] = [];
  @Input() clusterCentralityAvg: ClusterCentralityAvg = {
    degree_centrality: 0,
    average_clustering: 0,
    closeness_centrality: 0,
  };
  /**
   * Pre-existing behavior being preserved, not introduced here: this always uses the *main*
   * view's ordinal color scale, even when rendering the compare instance — the same systemic
   * "categorical fills always use colorScale, never colorScaleCompare" bug already documented
   * on HexagonRenderContext.colorScale in HexagonViewComponent. Confirmed via the original
   * template (both main and compare blocks call `colorScale(...)`, never `colorScaleCompare`).
   */
  @Input() colorScale!: (value: string) => string;

  /**
   * Renders the Plotly neighborhood-enrichment bar chart for `leiden`'s cluster into this
   * instance's own `#cluster-nhood-heatmap`/`#cluster-nhood-heatmap-compare` div. Called by the
   * parent's onTabChange (still parent-resident — see class doc) via @ViewChild, with explicit
   * parameters rather than @Input()s, since it's invoked imperatively right after the parent
   * mutates its own state, in the same tick Angular's own change-detection hasn't run yet.
   */
  public renderNhoodHeatmap(leiden: number | null | undefined, leidenClusterAnnotations: Record<string, any> | undefined): void {
    const clusterAnnotation = leiden === null || leiden === undefined ? null : leidenClusterAnnotations?.[String(leiden)];
    const enrichment = clusterAnnotation?.neighborhood_enrichment;
    const containerId = this.isCompare ? 'cluster-nhood-heatmap-compare' : 'cluster-nhood-heatmap';

    if (!enrichment || !Array.isArray(enrichment)) return;

    const n = enrichment.length;
    const clusterLabels = Array.from({ length: n }, (_, i) => `Cluster ${i}`);

    const minValue = Math.min(...enrichment);
    const maxValue = Math.max(...enrichment);
    const normalized = (maxValue > minValue)
      ? enrichment.map((v: number) => (v - minValue) / (maxValue - minValue))
      : enrichment.map(() => 0);

    const data: Partial<Plotly.PlotData>[] = [
      {
        x: clusterLabels,
        y: normalized,
        type: 'bar',
        marker: { color: 'rgba(55, 128, 191, 0.7)' },
        name: `Cluster ${leiden} Neighborhood Enrichment`,
      }
    ];

    const layout = {
      margin: { t: 30, l: 60, r: 10, b: 40 },
      width: 300,
      height: 170,
      xaxis: {
        title: { text: 'Cluster' },
        automargin: true,
        tickfont: { size: 10 },
      },
      yaxis: {
        title: { text: 'Enrichment' },
        automargin: true,
        tickfont: { size: 10 },
      },
    };

    const container = document.getElementById(containerId);
    if (!container) {
      console.error(`Container ${containerId} not found for rendering heatmap`);
      return;
    }

    Plotly.purge(container);
    Plotly.newPlot(container, data, layout, { displayModeBar: false });
  }
}
