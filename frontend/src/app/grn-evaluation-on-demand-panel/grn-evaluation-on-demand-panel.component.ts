import { Component, EventEmitter, Input, Output } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FilterableTableComponent } from '../filterable-table/filterable-table.component';
import { CellFeature } from '../hexagon-view/cell-feature.types';
import { MatExpansionModule } from '@angular/material/expansion';
import { MatFormField, MatLabel } from '@angular/material/form-field';
import { MatOption } from '@angular/material/autocomplete';
import { MatSelect } from '@angular/material/select';
import { MatInputModule } from '@angular/material/input';
import { MatCheckboxModule } from '@angular/material/checkbox';
import { MatButtonModule } from '@angular/material/button';
import { MatIconModule } from '@angular/material/icon';
import { MatProgressSpinnerModule } from '@angular/material/progress-spinner';

type GrnTableData = { [col: string]: { [index: string]: string | number } } | string[];
/** Matches GrnEvaluationPanelComponent's own GrnGraphFilterMode — duplicated rather than imported
 * so the two sibling panel components stay independently reusable, consistent with how the other
 * extracted sidebar-tab panels in this refactor (e.g. DgeaPanelComponent/FootprintPanelComponent)
 * don't cross-import from one another even where the parent happens to share state between them. */
export type GrnGraphFilterMode = 'prior' | 'extended' | 'full';

/**
 * Renders the "GRN Evaluation - On Demand" sidebar tab for one side (main or compare) — the
 * "Compute New GRN Evaluation" form, the "Load Existing GRN Evaluation Results" run picker, and
 * (once a run is selected) its Statistics/Graph/Plots results plus a "Register Run to Dataset"
 * button. Note the tab label is the identical literal string ("GRN Evaluation - On Demand") on
 * both the main and compare `mat-tab-group`s — the two instances are told apart by which
 * tab-group they are nested under, not by their label text.
 *
 * A genuinely passive, @Input()-driven component, for the same reasons as
 * GrnEvaluationPanelComponent (see its class doc): every fetch (computeOnDemandGrnEvaluation,
 * loadOnDemandGrnRuns, onSelectGrnRun => loadGrnGraph/loadGrnStats/loadGrnPlots,
 * registerOnDemandGrnRun) and D3-render (visualizeGrnGraph — looked up by the parent via a global
 * `document.querySelector`/`d3.select`, not through this component) stays in
 * HexagonPlotComponent.
 *
 * Bugs found while extracting this component (all found, NOT fixed — fixing any of them requires
 * editing HexagonPlotComponent's fields/methods, which is out of scope for this extraction):
 *
 * 1. The "Compute New GRN Evaluation" form's inputs (Observation Key/Cluster/Name/Gene Set/Use
 *    Prior GRN) and its error message are the ONLY pieces of this tab's state that are NOT split
 *    into main/compare pairs: `onDemandGrnParams`/`onDemandGrnError` are single fields shared by
 *    both the main and compare "Compute New GRN Evaluation" forms — unlike `onDemandOutputPath`/
 *    `onDemandOutputPathCompare`, `onDemandGrnRuns`/`onDemandGrnRunsCompare`,
 *    `selectedOnDemandRun`/`selectedOnDemandRunCompare`, `grnTabError`/`grnTabErrorCompare`, and
 *    every other piece of this tab's state, which ARE all correctly split (confirmed by grepping
 *    the whole file for both names — there is no `onDemandGrnParamsCompare`/
 *    `onDemandGrnErrorCompare` declared anywhere). So today, typing into either side's compute
 *    form fills in the SAME obsKey/cluster/name/geneSet/usePriorGrn values on both sides, and an
 *    error from either side's compute request shows on both. This is the same shape of bug this
 *    refactor has repeatedly found elsewhere (a `compare` parameter accepted but the underlying
 *    data not actually kept separate): computeOnDemandGrnEvaluation(compare) DOES correctly use
 *    `compare` for the dataset lookup and the output-path field, but reads
 *    `this.onDemandGrnParams`/writes `this.onDemandGrnError` unconditionally regardless of
 *    `compare`. This component's `obsKey`/`cluster`/`name`/`geneSet`/`usePriorGrn`/`computeError`
 *    @Input()s are independent, generic per-instance values (the usual convention for this
 *    refactor) — the parent CAN fix this by introducing genuinely separate
 *    `onDemandGrnParamsCompare`/`onDemandGrnErrorCompare` fields and binding them to the compare
 *    instance, but wiring both instances to today's same shared fields will preserve the bug
 *    exactly as it is today.
 * 2. Informational only, not a main/compare bug (symmetric on both sides) — the "GRN Graph"
 *    results section is entirely non-functional: its container div
 *    (`#grn_gene_net_container(_compare)`) is never targeted by any rendering code.
 *    onSelectGrnRun's loadGrnGraph fetches all 3 filter modes and, via
 *    visualizeGrnGraphForFilter, calls visualizeGrnGraph(graphData, compare) with no explicit
 *    filterMode argument, which makes visualizeGrnGraph fall back to reading
 *    `this.grnGraphFilter(Compare)` and target `#tf_gene_net_container_<mode>(_compare)` instead —
 *    i.e. the *pre-computed* "GRN Evaluation" tab's own containers (see
 *    GrnEvaluationPanelComponent's bug (2)). Relatedly, this section's
 *    `isLoadingGrnGraph`/`isLoadingGrnGraphCompare` flag (gating its spinner and the container's
 *    `display: none`/`block`) is initialized to `false` and never reassigned anywhere in the
 *    file, so the spinner never shows and the (permanently empty) container is always visible.
 *    This component's `isLoadingGraph`/graph container faithfully reproduce this dead state
 *    rather than inventing a fix for it.
 * 3. Informational only, not fixed (shared with GrnEvaluationPanelComponent's bug (3)): this tab's
 *    Graph Filter `<mat-select>` calls the same onChangeGrnGraphFilter(value, compare) as the
 *    pre-computed tab's filter selector, whose "no cached data yet" branch is a no-op — so
 *    switching filters here before loadGrnGraph's fetch for that mode resolves silently shows
 *    nothing new.
 */
@Component({
  selector: 'app-grn-evaluation-on-demand-panel',
  standalone: true,
  imports: [
    CommonModule,
    FilterableTableComponent,
    MatExpansionModule,
    MatFormField,
    MatLabel,
    MatOption,
    MatSelect,
    MatInputModule,
    MatCheckboxModule,
    MatButtonModule,
    MatIconModule,
    MatProgressSpinnerModule,
  ],
  templateUrl: './grn-evaluation-on-demand-panel.component.html',
  styleUrls: ['./grn-evaluation-on-demand-panel.component.scss'],
})
export class GrnEvaluationOnDemandPanelComponent {
  @Input() isCompare = false;

  // --- Compute New GRN Evaluation form (see bug (1) above) ---
  @Input() obsKey = 'cell_type';
  @Output() obsKeyChange = new EventEmitter<string>();
  @Input() cluster = '0';
  @Output() clusterChange = new EventEmitter<string>();
  @Input() name = 'GRN_Evaluation';
  @Output() nameChange = new EventEmitter<string>();
  @Input() geneSet = '';
  @Output() geneSetChange = new EventEmitter<string>();
  @Input() usePriorGrn = false;
  @Output() usePriorGrnChange = new EventEmitter<boolean>();

  @Input() isComputingOnDemandGrn = false;
  /** Parent is expected to call its own computeOnDemandGrnEvaluation(isCompare) in response. */
  @Output() computeRequested = new EventEmitter<void>();
  @Input() computeError = '';
  @Input() outputPath: string | null = null;

  get computeFormTitle(): string {
    return this.isCompare ? 'Compute New GRN Evaluation (Compare)' : 'Compute New GRN Evaluation';
  }

  // --- Load Existing GRN Evaluation Results ---
  @Input() availableRuns: string[] = [];
  @Input() selectedRun: string | null = null;
  /** Parent is expected to call its own onSelectGrnRun(value, isCompare) in response — not just
   * record the new value — since that also triggers loadGrnGraph/loadGrnStats/loadGrnPlots. */
  @Output() selectedRunChange = new EventEmitter<string>();
  /** Parent is expected to call its own loadOnDemandGrnRuns(isCompare) in response. */
  @Output() refreshRunsRequested = new EventEmitter<void>();
  @Input() tabError: string | null = null;

  // --- Results: Statistics ---
  @Input() isLoadingStats = false;
  @Input() peakStats: GrnTableData | null = null;
  @Input() motifStats: GrnTableData | null = null;
  @Input() features: CellFeature[] = [];
  @Input() datasetId?: string;
  @Output() tableLoadingChange = new EventEmitter<{ loading: boolean; isCompare: boolean }>();

  // --- Results: GRN Graph (see bug (2) above — this section is currently non-functional) ---
  @Input() graphFilter: GrnGraphFilterMode = 'full';
  /** Parent is expected to call its own onChangeGrnGraphFilter(value, isCompare) in response. */
  @Output() graphFilterChange = new EventEmitter<GrnGraphFilterMode>();
  @Input() isLoadingGraph = false;

  get graphContainerId(): string {
    return this.isCompare ? 'grn_gene_net_container_compare' : 'grn_gene_net_container';
  }

  // --- Results: Plots ---
  @Input() isLoadingPlots = false;
  @Input() plots: { [directory: string]: string[] } = {};
  /** Parent is expected to bind `(name) => this.getGrnImageUrl(name, isCompare)`. */
  @Input() getImageUrl: (imageName: string) => string = () => '';

  /** Parent is expected to call its own registerOnDemandGrnRun(isCompare) in response. */
  @Output() registerRequested = new EventEmitter<void>();
}
