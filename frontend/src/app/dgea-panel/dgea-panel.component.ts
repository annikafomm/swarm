import { Component, EventEmitter, Input, Output } from '@angular/core';
import { CommonModule } from '@angular/common';
import * as Plotly from 'plotly.js-dist-min';
import { FilterableTableComponent } from '../filterable-table/filterable-table.component';
import { CellFeature } from '../hexagon-view/cell-feature.types';
import { MatFormField, MatLabel } from '@angular/material/form-field';
import { MatOption } from '@angular/material/autocomplete';
import { MatSelect } from '@angular/material/select';
import { MatCheckboxModule } from '@angular/material/checkbox';

type DgeaTableData = { [col: string]: { [index: string]: string | number } } | string[];

export interface DgeaHeatmapContext {
  groups: (string | number)[];
  rows: { gene: string; scaled: number; raw: number }[];
}

export interface DgeaComparison {
  group1: string;
  group2: string;
  n1: number;
  n2: number;
  skipped?: boolean;
  skip_reason?: string;
  table: DgeaTableData;
  heatmap_context?: DgeaHeatmapContext | null;
}

/**
 * Renders the "DGEA"/"Compare - DGEA" sidebar tab for one side (main or compare) — the
 * Grouping/Group1/Group2/vs-All selection UI, the resulting comparison's `<app-table>`, and the
 * context-heatmap container div. Extracted from HexagonPlotComponent, which used to
 * hand-duplicate this markup verbatim for main vs. compare.
 *
 * A genuinely passive, @Input()-driven component, like RegulatoryTablesPanelComponent: the actual
 * DGEA selection state (selectedDgeaObsCol/Group1/Group2/dgeaVsAll, per side), its mutation
 * handlers (onDgeaObsColChange/onDgeaSelectionChange/onDgeaVsAllChange/initDgeaSelection), and the
 * "show this gene on the map" cross-cutting behavior (onDgeaGeneSelected/showDgeaGeneOnMainPlot)
 * all stay in HexagonPlotComponent, since they also touch the color-by/gene-expression state that
 * belongs to later refactor phases. This component only emits `*Change`/`geneSelected` and the
 * parent's existing handlers absorb them — same shape as RegulatoryTablesPanelComponent's
 * `itemSelected`/`tableLoadingChange`.
 *
 * `renderDgeaHeatmap` mirrors ClusterInfoPanelComponent.renderNhoodHeatmap exactly: it's called
 * imperatively by the parent (via `@ViewChild('mainDgeaPanel'|'compareDgeaPanel')`) right after
 * the parent mutates its own DGEA selection state, in the same tick — not reactively through
 * @Input() change detection — so it takes explicit parameters for the two pieces of data the
 * original method read off `this` (the resolved heatmap context and the axis label), and looks up
 * its own container by id rather than by `@ViewChild`, for the same reason.
 *
 * Bugs found in the original markup while extracting this component:
 *
 * 1. FIXED here (self-contained to this component's own template): FilterableTableComponent's
 *    `show_on_plot` action — DGEA's only action column — emits via its `geneSelectedCompare`
 *    output instead of `geneSelected` whenever `isCompare` is true (see
 *    FilterableTableComponent.onShowAction), and never falls through to its own
 *    fetchAndUpdate/featuresUpdated (unlike every other action column), so there is no
 *    fallback path that still colors the map. The original "Compare - DGEA" `<app-table>` only
 *    listened on `(geneSelected)`, so clicking "show on plot" for a compare-side DGEA gene did
 *    nothing at all. Fixed by listening on both `(geneSelected)` and `(geneSelectedCompare)` here
 *    and forwarding both through this component's single `geneSelected` output — only one of the
 *    two ever fires for a given `isCompare`, so this is safe for the main instance too.
 * 2. FIXED here (structural side effect, not a behavior change): the original main and compare
 *    heatmap divs both used the literal DOM id `dgea-heatmap` (the compare div's `id` should have
 *    read `dgea-heatmap-compare`, matching its own `#dgeaHeatmapCompare` template-ref/`@ViewChild`
 *    and Plotly target). This component derives its container id from `isCompare`
 *    (`dgea-heatmap`/`dgea-heatmap-compare`), the same convention already used by
 *    ClusterInfoPanelComponent's `cluster-nhood-heatmap`/`cluster-nhood-heatmap-compare`, so the
 *    duplicate id can no longer occur.
 * 3. NOT fixed — would require editing HexagonPlotComponent's class/template, out of scope for
 *    this extraction (same situation as FootprintPanelComponent's documented
 *    footprintPlotUrls/availableMotifs cross-contamination bug):
 *    a. In the original "Compare - DGEA" tab, the Grouping/Group1/Group2/vs-All controls are
 *       `[(ngModel)]`-bound to the *main*-side `selectedDgeaObsCol`/`selectedDgeaGroup1`/
 *       `selectedDgeaGroup2`/`dgeaVsAll` fields, and its dropdown option lists call
 *       `getDgeaObsCols()`/`getDgeaLevels()`/`getSelectedDgeaComparison()` with no `compare`
 *       argument (defaulting to `false`, i.e. main's `meta`) — even though the *change handlers*
 *       on those same controls do pass `compare=true`. So today, using the Compare tab's dropdowns
 *       actually reads/mutates the *main* tab's selection state, while the dedicated
 *       `selectedDgeaObsColCompare`/`Group1Compare`/`Group2Compare`/`dgeaVsAllCompare` fields
 *       (already declared on HexagonPlotComponent for exactly this purpose) sit mostly unused.
 *    b. The original "Compare - DGEA" `<app-table>` also bound `[features]="features"` (main's
 *       array) and `[datasetId]="selectedDataset?.id"` (main's dataset) instead of
 *       `[features]="compareFeatures"` / `[datasetId]="selectedDatasetCompare?.id ||
 *       selectedDataset?.id"` — the pattern every other compare-side `<app-table>` in the file
 *       uses.
 *    This component's `@Input()`s (`selectedObsCol`, `group1`, `group2`, `vsAll`, `obsCols`,
 *    `levels`, `comparison`, `features`, `datasetId`) are generic, per-instance values — the same
 *    convention as CoOccurrencePanelComponent/FootprintPanelComponent — so whoever wires the
 *    compare instance CAN supply the correct compare-scoped values; simply wiring this component
 *    in without also passing them will preserve today's cross-contamination.
 */
@Component({
  selector: 'app-dgea-panel',
  standalone: true,
  imports: [CommonModule, FilterableTableComponent, MatFormField, MatLabel, MatOption, MatSelect, MatCheckboxModule],
  templateUrl: './dgea-panel.component.html',
  styleUrls: ['./dgea-panel.component.scss'],
})
export class DgeaPanelComponent {
  @Input() isCompare = false;

  /** dgeaReady / dgeaReadyCompare — whether this side has any DGEA data at all. */
  @Input() ready = false;

  /** getDgeaObsCols(isCompare) — the available grouping columns (e.g. 'cell_type', 'leiden'). */
  @Input() obsCols: string[] = [];
  /** selectedDgeaObsCol / selectedDgeaObsColCompare. */
  @Input() selectedObsCol = 'cell_type';
  /** Parent is expected to assign this back onto its own selectedDgeaObsCol(Compare) field and
   * call its existing onDgeaObsColChange(isCompare) — not just record the new value — since that
   * also resets group1/group2/vsAll and re-initializes the selection. */
  @Output() obsColChange = new EventEmitter<string>();

  /** getDgeaLevels(isCompare) — the available levels for Group 1/Group 2. */
  @Input() levels: string[] = [];
  /** selectedDgeaGroup1 / selectedDgeaGroup1Compare. */
  @Input() group1: string | null = null;
  /** Parent absorbs via its onDgeaSelectionChange(isCompare), same as group2Change. */
  @Output() group1Change = new EventEmitter<string>();
  /** selectedDgeaGroup2 / selectedDgeaGroup2Compare. */
  @Input() group2: string | null = null;
  @Output() group2Change = new EventEmitter<string>();
  /** dgeaVsAll / dgeaVsAllCompare. */
  @Input() vsAll = false;
  /** Parent absorbs via its onDgeaVsAllChange(isCompare). */
  @Output() vsAllChange = new EventEmitter<boolean>();

  /** getSelectedDgeaComparison(isCompare) — null renders the "no comparison found" empty state. */
  @Input() comparison: DgeaComparison | null = null;

  @Input() features: CellFeature[] = [];
  @Input() datasetId?: string;

  @Output() tableLoadingChange = new EventEmitter<{ loading: boolean; isCompare: boolean }>();
  /**
   * Emits for the DGEA table's only action, "show on plot" — see bug (1) above for why this
   * listens on both of `<app-table>`'s `geneSelected`/`geneSelectedCompare` outputs. The parent is
   * expected to call its own onDgeaGeneSelected(event, isCompare) in response.
   */
  @Output() geneSelected = new EventEmitter<{ gene: string; action: string }>();

  getObsColLabel(col: string): string {
    if (col === 'cell_type') return 'Cell type';
    if (col === 'leiden') return 'Leiden';
    return col;
  }

  get group1Label(): string {
    return this.selectedObsCol === 'cell_type' ? 'Cell type 1' : 'Group 1';
  }

  get group2Label(): string {
    return this.selectedObsCol === 'cell_type' ? 'Cell type 2' : 'Group 2';
  }

  /**
   * Renders the DGEA context heatmap into this instance's own `#dgea-heatmap`/
   * `#dgea-heatmap-compare` div. Called by the parent (its onTabChange, and after
   * onDgeaObsColChange/onDgeaSelectionChange/onDgeaVsAllChange) via `@ViewChild`, with explicit
   * parameters — the heatmap context the parent's getSelectedDgeaHeatmap(isCompare) resolved, and
   * the axis label its getDgeaObsColLabel(...) resolved — rather than @Input()s, since it runs
   * before Angular's own change detection would have propagated those as bindings.
   */
  public renderDgeaHeatmap(heatmap: DgeaHeatmapContext | null, xAxisTitle: string): void {
    const containerId = this.isCompare ? 'dgea-heatmap-compare' : 'dgea-heatmap';
    const container = document.getElementById(containerId);
    if (!container) return;

    if (!heatmap || !heatmap.groups || !heatmap.rows || !heatmap.rows.length) {
      Plotly.purge(container);
      return;
    }

    const x = heatmap.groups.map(g => String(g));
    const y = heatmap.rows.map(r => r.gene);
    const z = heatmap.rows.map(r => r.scaled);
    const raw = heatmap.rows.map(r => r.raw);

    const data: Partial<Plotly.PlotData>[] = [
      {
        type: 'heatmap',
        x,
        y,
        z,
        customdata: raw,
        colorscale: 'RdBu',
        reversescale: true,
        hovertemplate:
          'Gene: %{y}<br>' +
          'Group: %{x}<br>' +
          'Scaled expression: %{z:.2f}<br>' +
          'Mean expression: %{customdata:.2f}<extra></extra>',
      }
    ];

    const layout: Partial<Plotly.Layout> = {
      margin: { t: 30, l: 140, r: 20, b: 100 },
      height: Math.max(420, y.length * 22),
      xaxis: {
        title: { text: xAxisTitle },
        tickangle: -45,
        automargin: true,
        type: 'category'
      },
      yaxis: {
        title: { text: 'Genes' },
        automargin: true,
        autorange: 'reversed'
      }
    };

    Plotly.purge(container);
    Plotly.newPlot(container, data, layout, {
      responsive: true,
      displayModeBar: false
    });
  }
}
