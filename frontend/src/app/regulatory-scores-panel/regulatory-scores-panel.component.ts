import { Component, EventEmitter, Input, Output } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FilterableTableComponent } from '../filterable-table/filterable-table.component';
import { CellFeature } from '../hexagon-view/cell-feature.types';
import { MatFormField, MatLabel } from '@angular/material/form-field';
import { MatOption } from '@angular/material/autocomplete';
import { MatSelect } from '@angular/material/select';
import { MatExpansionModule } from '@angular/material/expansion';
import { MatProgressSpinnerModule } from '@angular/material/progress-spinner';
import { MatIconModule } from '@angular/material/icon';
import { MatButtonModule } from '@angular/material/button';
import { MatTooltipModule } from '@angular/material/tooltip';
import { InfoService } from '../info.service';

type TableData = { [col: string]: { [index: string]: string | number } } | string[];

export interface RegulatorySliderData {
  step: number;
  min_border: number;
  max_border: number;
  default_value: number;
}

/**
 * Renders the "Regulatory Scores"/"Compare - Regulatory Scores" sidebar tab for one side (main or
 * compare) — the Genie3/SpongEffects `<app-table>`s (both the dataset-wide "Global Regulatory
 * Scores" pair and the per-cell "Cell Regulatory Scores" pair), the Genie3/Sponge gene-set
 * dropdowns, their weight/p-value/min-edges sliders, and the plain D3-target container divs for
 * the two network graphs. Extracted from HexagonPlotComponent, which used to hand-duplicate this
 * markup verbatim for main vs. compare. This does NOT cover "GRN Evaluation" (peak_stats/
 * motif_stats/tf_gene_net_container), a different tab extracted separately.
 *
 * What stays in HexagonPlotComponent, and why (every call site checked):
 * - `updateSubgraphGenie3`/`updateSubgraphSponge` (the actual HTTP fetch of the network + slider
 *   bounds, using the same `requestTokens`/`nextRequestToken` stale-request guard already relied on
 *   elsewhere — left untouched) and `visualizeGenie3Subgraph`/`visualizeSpongeSubgraph` (the D3
 *   force-graph render into `#aucell_graph_genie3(_compare|_sponge...)`) are called from far more
 *   call sites than just this tab's own controls: `onGeneSetChange`, cell clicks
 *   (`displayCellDetails`), cluster selection (`selectCluster`, `selectCellFromCluster`,
 *   `autoSelectDefaultCluster`), and the panel-open hook (`onRegulatoryGraphsPanelOpened`). Per the
 *   established convention (DgeaPanelComponent/FootprintPanelComponent), multi-call-site
 *   fetch/render logic stays in the parent; this component only supplies the container div (by the
 *   exact same `id` the parent's `d3.select(...)` already targets) and emits request events.
 * - `getRegulatoryTableData`/`getRegulatoryActionColumns`/`getRegulatoryEmptyMessage` are resolved
 *   by the parent (also driven by `regulatoryScoreDisplayMode(Compare)`, itself changed from
 *   multiple places) and passed in as plain `@Input()`s — this component does not know about
 *   `meta`/`metaCompare` at all.
 * - `onGeneSelectedFromTable`/`selectGeneSetFromTable`/`onGeneSetChange` mutate
 *   `selectedRegulatoryScore(Compare)`/`selectedGeneSet*(Compare)`/`previousGeneSet*(Compare)` and
 *   then trigger the fetch/render above — kept in the parent since they're the same functions the
 *   cross-cutting call sites above also invoke.
 * - `getRegulatoryScoresforSpots` (called from `displayCellDetails`) populates
 *   `genie3RawData(Compare)`/`spongeRawData(Compare)`/`genie3Elements(Compare)`/
 *   `spongeElements(Compare)` — also left in the parent as a multi-purpose, multi-call-site fetch.
 *
 * Bugs found in the original markup and FIXED here (safe: confined to markup this component now
 * owns, no parent-file change required):
 * 1. The Compare tab's Genie3/Sponge "Gene Set" `<mat-select>`s iterated over `genie3Elements`/
 *    `spongeElements` (the MAIN side's element lists) instead of `genie3ElementsCompare`/
 *    `spongeElementsCompare` — so the Compare network-graph gene-set dropdowns actually offered the
 *    Main dataset's regulators. Fixed by making `genie3Elements`/`spongeElements` plain per-instance
 *    `@Input()`s (same convention as every sibling panel); whoever wires the compare instance must
 *    supply `genie3ElementsCompare`/`spongeElementsCompare` — simply wiring this component in without
 *    doing so would silently reproduce today's bug.
 * 2. Main's two "graph-summary" lines (`*ngIf="!isLoadingGenie3 && genie3WeightCutoff !== null"` /
 *    `*ngIf="!isLoadingSponge && spongePValueCutoff !== null"`) used a condition that is always true
 *    once not loading, since `genie3WeightCutoff`/`spongePValueCutoff` are plain numbers that are
 *    never `null` — so the "Genie3 Weight Cutoff: … · Number of edges: …" summary rendered even with
 *    no gene set selected, right alongside the "No Genie3 network selected" message. The Compare
 *    side already used the correct condition (`selectedGeneSetGenie3Compare`/
 *    `selectedGeneSetSpongeCompare`). Fixed by using `selectedGeneSetGenie3`/`selectedGeneSetSponge`
 *    uniformly, matching Compare's already-correct behavior.
 * 3. Compare's "Global Regulatory Scores" and "Regulatory Network Graphs" `<mat-expansion-panel>`s
 *    were not wrapped in a `<mat-accordion>` at all (Main wraps all three of its panels in one).
 *    Fixed by wrapping all three panels (Global / Graphs / Cell) in one shared `<mat-accordion>` for
 *    both instances — the only structural option once one physical template serves both sides.
 * 4. Compare's two "Cell Regulatory Scores" `<app-table>`s (Genie3 and Sponge) omitted
 *    `[emptyMessage]` entirely, unlike Main's `'No Cell Regulatory Scores Loaded'` — the same class
 *    of compare-side omission RegulatoryTablesPanelComponent already documented and fixed for its
 *    four tabs. Fixed by supplying it uniformly.
 * Minor cosmetic normalizations (not behavioral bugs, just picks made necessary by unifying one
 * template): Main's per-cell Sponge heading was `<h4>` while its Genie3 sibling and both of
 * Compare's headings were `<h3>` — normalized to `<h3>`. The Global Scores panel's `expanded="true"`
 * (Main, plain attribute) vs `[expanded]="true"` (Compare, bound) is normalized to `[expanded]="true"`.
 *
 * Bugs found, NOT fixed (live inside methods that stay in HexagonPlotComponent, or require changing
 * its class/template — out of scope for this file-only extraction; flagged in detail for whoever
 * wires this component in):
 * a. `visualizeGenie3Subgraph(compare)` slices edges with `candidateEdges.slice(0, this.genie3MinEdges)`
 *    — ignores `compare` entirely, so the Compare graph always uses Main's `genie3MinEdges`, never
 *    `genie3MinEdgesCompare`.
 * b. `visualizeSpongeSubgraph(compare)`, inside its own `if (compare)` branch, filters edges with
 *    `edge.p_adjusted < this.spongePValueCutoff` — reads Main's cutoff instead of
 *    `spongePValueCutoffCompare` even though the surrounding branch is specifically for compare.
 * c. `visualizeSpongeSubgraph(compare)` also slices with `candidateEdges.slice(0, this.spongeMinEdges)`
 *    — same unconditional-Main-value bug as (a), for Sponge.
 * d. `visualizeSpongeSubgraph` computes `targets` from `this.geneSetsGenie3[regulator]` instead of
 *    `this.geneSetsSponge[regulator]` — a Genie3/Sponge (not main/compare) copy-paste bug that breaks
 *    node-coloring (regulator's targets almost never actually match) for both main and compare Sponge
 *    graphs equally.
 * e. `visualizeGenie3Subgraph`/`visualizeSpongeSubgraph` read `this.genie3Width`/`this.spongeWidth`
 *    unconditionally (no `compare` branch at all) to size the SVG, and there is no
 *    `genie3WidthCompare`/`spongeWidthCompare`. `updateGraphWidths()` only measures
 *    `@ViewChild('aucell_graph_genie3'|'aucell_graph_sponge')` — the Main-only template refs — so the
 *    Compare graph's width is always derived from Main's container (or the 600px default), never its
 *    own. Relatedly: once this div's markup lives in this child component, those two `@ViewChild`
 *    lookups will return `undefined` for good (a parent's `@ViewChild` cannot reach into a child
 *    component's own template) — `updateGraphWidths()` will need reworking regardless (e.g. reading
 *    `document.getElementById(...)` by the same `id`s this component still renders) to keep working
 *    at all, not just to fix the compare-width bug.
 * f. `analyzeGeneSetInGProfiler()` ignores which button was clicked and is hardcoded to
 *    `this.selectedGeneSetGenie3`/`this.geneSetsGenie3` — so "Analyze All Sponge Genes in gProfiler"
 *    actually analyzes the Genie3 gene set. It also has no `compare` parameter and today's Compare tab
 *    has no gProfiler button at all. This component faithfully reproduces that (the button only
 *    renders when `!isCompare`, via `analyzeInGProfilerRequested`), but the emitted event does carry
 *    the network type so the wiring step can fix the mixup if it also updates
 *    `analyzeGeneSetInGProfiler` to accept it.
 * g. `onRegulatoryGraphsPanelOpened()` is hardcoded main-only (no `compare` parameter) and is wired
 *    only via `(opened)` on Main's "Regulatory Network Graphs" panel; Compare's panel has no `(opened)`
 *    handler at all, so — combined with `matExpansionPanelContent`'s lazy content — a pre-selected
 *    Compare gene set's graph may not actually render until the panel is opened AND the gene set is
 *    re-selected. This component emits `graphsPanelOpened` from both instances (harmless by itself —
 *    an unbound `@Output()` is a no-op); binding it on the compare instance too, to fix the gap, first
 *    requires giving `onRegulatoryGraphsPanelOpened` a `compare` parameter in the parent.
 * h. Compare's sliders are wired to `(input)="onGeneSetChange(true)"` instead of a dedicated
 *    slider-changed handler. `onGeneSetChange` only re-renders when the *gene set* value itself
 *    changed, not the slider values, so dragging Compare's weight/min-edges/p-value sliders today
 *    triggers no re-render at all (only changing the Gene Set dropdown does). Main's sliders correctly
 *    call `onSliderChangeGenie3()`/`onSliderChangeSponge()` — but those two are themselves hardcoded
 *    main-only (no `compare` parameter), so simply binding this component's new
 *    `genie3SliderChanged`/`spongeSliderChanged` outputs to them for the compare instance would
 *    incorrectly refresh Main's graph instead. Needs a `compare`-aware slider handler in the parent.
 * i. `requestTokens` (used by `nextRequestToken`/`updateSubgraphGenie3`/`updateSubgraphSponge`) is
 *    keyed only by `'genie3'`/`'sponge'`, not split by side — a fresh Main request and a fresh Compare
 *    request for the same network type share one stale-request counter, so one side's new request can
 *    mark the other side's in-flight request "stale" and suppress its loading-flag reset. Symmetric
 *    between main/compare (not a copy-paste asymmetry) and explicitly out of scope to touch here.
 * j. `genie3Elements`/`spongeElements` (and their Compare counterparts) do double duty: populated once
 *    from `geneSetsGenie3`/`geneSetsSponge` keys at dataset load (for the Gene Set dropdown — see
 *    fixed bug 1 above), then silently overwritten by `getRegulatoryScoresforSpots` on every cell click
 *    with a different set of names (the per-cell score-column keys). After clicking any cell, the Gene
 *    Set dropdown's options change to the wrong list. Pre-existing, symmetric between main/compare, out
 *    of scope.
 * k. `globalRegulatoryScoresTutorial()`'s Shepherd step targets `#global-regulatory-scores-info`, an
 *    element id that does not exist anywhere in the current template (the button only carries class
 *    `info-btn`) — the tour never actually anchors to the button. Not a main/compare bug, so left
 *    unfixed; if ever fixed, the `id` would need to be added inside this component's own template,
 *    since the button now lives here.
 */
@Component({
  selector: 'app-regulatory-scores-panel',
  standalone: true,
  imports: [
    CommonModule,
    FilterableTableComponent,
    MatFormField,
    MatLabel,
    MatOption,
    MatSelect,
    MatExpansionModule,
    MatProgressSpinnerModule,
    MatIconModule,
    MatButtonModule,
    MatTooltipModule,
  ],
  templateUrl: './regulatory-scores-panel.component.html',
  styleUrls: ['./regulatory-scores-panel.component.scss'],
})
export class RegulatoryScoresPanelComponent {
  constructor(public infoService: InfoService) {}

  @Input() isCompare = false;

  // ---- Current-selection summary (also reused by the graphs section below) ----
  /** selectedRegulatoryScore / selectedRegulatoryScoreCompare. */
  @Input() selectedRegulatoryScore: string | null = null;

  // ---- Global Regulatory Scores ----
  /** regulatoryScoreDisplayMode / regulatoryScoreDisplayModeCompare. */
  @Input() regulatoryScoreDisplayMode: 'raw' | 'moranI' | 'gearyC' = 'raw';
  /** Parent is expected to assign this back onto its own field and call its existing
   * onRegulatoryDisplayModeChange(isCompare). */
  @Output() regulatoryScoreDisplayModeChange = new EventEmitter<'raw' | 'moranI' | 'gearyC'>();

  /** getRegulatoryTableData('genie3'|'sponge', isCompare). */
  @Input() genie3TableData: TableData = {};
  @Input() spongeTableData: TableData = {};
  /** getRegulatoryActionColumns('genie3'|'sponge', isCompare). */
  @Input() genie3ActionColumns: string[] = [];
  @Input() spongeActionColumns: string[] = [];
  /** getRegulatoryEmptyMessage('genie3'|'sponge', isCompare). */
  @Input() genie3EmptyMessage = '';
  @Input() spongeEmptyMessage = '';

  @Input() features: CellFeature[] = [];
  @Input() datasetId?: string;

  @Output() tableLoadingChange = new EventEmitter<{ loading: boolean; isCompare: boolean }>();
  /**
   * Shared by all four `<app-table>`s in this component (both "Global" tables and both "Cell"
   * tables). Listens on both `geneSelected`/`geneSelectedCompare` — only one ever fires for a given
   * `isCompare` — and adds `networkType` so the parent's existing
   * onGeneSelectedFromTable(event, networkType, isCompare) can be called unmodified.
   */
  @Output() geneSelected = new EventEmitter<{ gene: string; action: string; networkType: 'genie3' | 'sponge' }>();
  /** Relayed from the two "Global Regulatory Scores" tables only. The parent binds
   * updateHexColors() (main) / repaintBothViews() (compare) — matching the original's pre-existing
   * main/compare asymmetry (see cellFeaturesUpdated below, which both sides bind to
   * repaintBothViews()). */
  @Output() globalFeaturesUpdated = new EventEmitter<void>();
  /** Relayed from the two "Cell Regulatory Scores" tables only. Both main and compare originally
   * bind this to repaintBothViews(), not updateHexColors(). */
  @Output() cellFeaturesUpdated = new EventEmitter<void>();

  /** Emits when the (main-only) help icon next to "Global Regulatory Scores" is clicked. Parent
   * calls its existing globalRegulatoryScoresTutorial(). */
  @Output() globalScoresTutorialRequested = new EventEmitter<void>();

  // ---- Regulatory Network Graphs: shared ----
  /** Emitted when the "Regulatory Network Graphs" panel is opened, from both instances. Binding
   * this on the compare instance to fix bug (g) above requires first giving
   * onRegulatoryGraphsPanelOpened a `compare` parameter in the parent — see class doc. */
  @Output() graphsPanelOpened = new EventEmitter<void>();

  // ---- Regulatory Network Graphs: Genie3 ----
  /** genie3Elements / genie3ElementsCompare — see FIXED bug 1 in the class doc: the compare
   * instance MUST be given genie3ElementsCompare here, not genie3Elements. */
  @Input() genie3Elements: string[] = [];
  /** selectedGeneSetGenie3 / selectedGeneSetGenie3Compare. */
  @Input() selectedGeneSetGenie3: string | null = null;
  /** Parent is expected to assign this back onto its own field and call its existing
   * onGeneSetChange(isCompare). */
  @Output() selectedGeneSetGenie3Change = new EventEmitter<string | null>();
  @Input() genie3SliderData: RegulatorySliderData | null = null;
  @Input() genie3WeightCutoff = 0.5;
  @Output() genie3WeightCutoffChange = new EventEmitter<number>();
  @Input() genie3MinEdges = 25;
  @Output() genie3MinEdgesChange = new EventEmitter<number>();
  /** Fired on every weight-cutoff/min-edges slider (input), after the *Change outputs above. See
   * NOT-FIXED bug (h): binding this correctly for the compare instance requires a compare-aware
   * slider handler in the parent, which doesn't exist today. */
  @Output() genie3SliderChanged = new EventEmitter<void>();
  @Input() isLoadingGenie3 = false;

  // ---- Regulatory Network Graphs: Sponge ----
  /** spongeElements / spongeElementsCompare — same FIXED-bug-1 caveat as genie3Elements. */
  @Input() spongeElements: string[] = [];
  /** selectedGeneSetSponge / selectedGeneSetSpongeCompare. */
  @Input() selectedGeneSetSponge: string | null = null;
  @Output() selectedGeneSetSpongeChange = new EventEmitter<string | null>();
  @Input() spongeSliderData: RegulatorySliderData | null = null;
  @Input() spongePValueCutoff = 0.05;
  @Output() spongePValueCutoffChange = new EventEmitter<number>();
  @Input() spongeMinEdges = 25;
  @Output() spongeMinEdgesChange = new EventEmitter<number>();
  @Output() spongeSliderChanged = new EventEmitter<void>();
  @Input() isLoadingSponge = false;

  /** Emits which network type's gProfiler button was clicked. Only ever fires from the main
   * instance — the button is hidden entirely when isCompare (see class doc bug f) — but carries
   * networkType regardless, in case the parent's analyzeGeneSetInGProfiler is ever fixed to use it. */
  @Output() analyzeInGProfilerRequested = new EventEmitter<'genie3' | 'sponge'>();

  // ---- Cell Regulatory Scores ----
  /** selectedCell / selectedCellCompare — also gates whether this panel renders at all. */
  @Input() selectedCell: CellFeature | null = null;
  /** genie3RawData / genie3RawDataCompare. */
  @Input() genie3RawData: TableData = {};
  /** spongeRawData / spongeRawDataCompare. */
  @Input() spongeRawData: TableData = {};

  readonly genie3CellActionColumns = [
    'aucell_scores_genie3',
    'viper_scores_genie3',
    'spongeffects_ssGSEA_scores_genie3',
    'spongeffects_GSVA_scores_genie3',
  ];
  readonly spongeCellActionColumns = [
    'aucell_scores_sponge',
    'viper_scores_sponge',
    'spongeffects_ssGSEA_scores_sponge',
    'spongeffects_GSVA_scores_sponge',
  ];

  get selectionLabel(): string {
    return this.isCompare ? 'Comparison selection' : 'Current selection';
  }

  /** #aucell_graph_genie3 / #aucell_graph_genie3_compare — the parent's D3 code still targets
   * these exact ids via `d3.select(...)`, which works regardless of which component renders the
   * div, since CSS ids are global in the DOM. */
  get genie3GraphContainerId(): string {
    return this.isCompare ? 'aucell_graph_genie3_compare' : 'aucell_graph_genie3';
  }
  get spongeGraphContainerId(): string {
    return this.isCompare ? 'aucell_graph_sponge_compare' : 'aucell_graph_sponge';
  }

  get minWeightSliderId(): string {
    return this.isCompare ? 'min-weight-slider-compare' : 'min-weight-slider';
  }
  get genie3MaxEdgesSliderId(): string {
    return this.isCompare ? 'max-edges-slider-compare' : 'max-edges-slider';
  }
  get maxPvalueSliderId(): string {
    return this.isCompare ? 'max-pvalue-slider-compare' : 'max-pvalue-slider';
  }
  get spongeMaxEdgesSliderId(): string {
    return this.isCompare ? 'max-edges-slider-sponge-compare' : 'max-edges-slider-sponge';
  }

  onGenie3WeightCutoffInput(value: string): void {
    this.genie3WeightCutoff = Number(value);
    this.genie3WeightCutoffChange.emit(this.genie3WeightCutoff);
    this.genie3SliderChanged.emit();
  }

  onGenie3MinEdgesInput(value: string): void {
    this.genie3MinEdges = Number(value);
    this.genie3MinEdgesChange.emit(this.genie3MinEdges);
    this.genie3SliderChanged.emit();
  }

  onSpongePValueCutoffInput(value: string): void {
    this.spongePValueCutoff = Number(value);
    this.spongePValueCutoffChange.emit(this.spongePValueCutoff);
    this.spongeSliderChanged.emit();
  }

  onSpongeMinEdgesInput(value: string): void {
    this.spongeMinEdges = Number(value);
    this.spongeMinEdgesChange.emit(this.spongeMinEdges);
    this.spongeSliderChanged.emit();
  }
}
