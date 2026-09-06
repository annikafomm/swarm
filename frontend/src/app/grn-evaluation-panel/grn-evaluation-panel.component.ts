import { Component, EventEmitter, Input, Output } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FilterableTableComponent } from '../filterable-table/filterable-table.component';
import { CellFeature } from '../hexagon-view/cell-feature.types';
import { MatExpansionModule } from '@angular/material/expansion';
import { MatFormField, MatLabel } from '@angular/material/form-field';
import { MatOption } from '@angular/material/autocomplete';
import { MatSelect } from '@angular/material/select';
import { MatProgressSpinnerModule } from '@angular/material/progress-spinner';

type GrnTableData = { [col: string]: { [index: string]: string | number } } | string[];
export type GrnGraphFilterMode = 'prior' | 'extended' | 'full';

/**
 * Renders the "GRN Evaluation"/"Compare - GRN Evaluation" sidebar tab for one side (main or
 * compare) — the pre-computed Peak-Gene Link Statistics / Motif Statistics `<app-table>`s, the
 * TF Gene Network Graph (filter-mode selector plus three D3-rendered container divs, one per
 * filter mode), and the GRN Evaluation Plots image gallery. Extracted from HexagonPlotComponent,
 * which used to hand-duplicate this markup verbatim for main vs. compare.
 *
 * A genuinely passive, @Input()-driven component — like FootprintPanelComponent, NOT like
 * DgeaPanelComponent/ClusterInfoPanelComponent's `renderX()`-via-@ViewChild pattern — because
 * every fetch/D3-render method behind this tab is invoked from more than one call site and/or
 * shares state with the sibling "GRN Evaluation - On Demand" tab (see
 * GrnEvaluationOnDemandPanelComponent's class doc), so per the refactor's "keep genuinely
 * cross-cutting orchestration in the parent" rule they all stay in HexagonPlotComponent:
 * - loadPrecomputedGrnGraph(compare) — fetches all 3 filter-mode graphs in parallel and stores
 *   them in grnGraphDataPrior/Extended/Full(Compare) — is called only from onTabChange's
 *   'GRN Evaluation'/'Compare - GRN Evaluation' branches (2 call sites total), but the D3 target
 *   container it visualizes into (via visualizeGrnGraph) is looked up by the *parent* with a
 *   plain `document.querySelector`/`d3.select('#id')`, not through this component — so no
 *   @ViewChild-exposed render method is needed here: this component only has to render the
 *   container divs with the right literal ids for the parent's existing code to find.
 * - onChangeGrnGraphFilter(newFilter, compare) — invoked from FOUR call sites: this tab's own
 *   filter <mat-select>, the sibling "GRN Evaluation - On Demand" tab's filter <mat-select>, and
 *   both of those again on the compare side — genuinely cross-cutting across two different
 *   sidebar tabs, not just main/compare.
 * - loadTfGraph(compare)/visualizeTfGraph — called from the same onTabChange branches as
 *   loadPrecomputedGrnGraph, but not wired to anything in the current template at all: it always
 *   targets `#tf_gene_net_container(_compare)`, an id that does not appear anywhere in
 *   hexagon-plot.component.html today (the three filter-mode containers below are
 *   `tf_gene_net_container_{prior,extended,full}(_compare)`). See "Bugs found" (2) below — this
 *   component exposes nothing for loadTfGraph, since its output has no visible target.
 * - loadPrecomputedGrnPlots(compare) — same "called only from onTabChange, 2 call sites" shape as
 *   loadPrecomputedGrnGraph, and shares its `grnPlots(Compare)`/`isLoadingGrnPlots(Compare)` state
 *   with the On-Demand tab's own Plots section (see bug (2) below).
 * - getGrnImageUrl(imageName, compare) — builds `<img>` srcs for both this tab's and the
 *   On-Demand tab's plot galleries; exposed here as a bound-function @Input(), the same
 *   convention ClusterInfoPanelComponent uses for `colorScale`.
 *
 * Bugs found while extracting this component (none are main/compare copy-paste bugs fixable
 * purely within this new component's own markup, so per the refactor's bug-handling convention
 * they are flagged here, not fixed):
 *
 * 1. FOUND, NOT FIXED HERE (fixable only by whoever wires this component into
 *    hexagon-plot.component.html, which is out of scope for this extraction): the *main*-side
 *    (non-compare) "Peak-Gene Link Statistics" and "Motif Statistics" `<app-table>`s currently
 *    bind `[features]="compareFeatures"` — the *compare*-side feature array — while every other
 *    main-side `<app-table>` in the file (and this tab's own compare instance) binds `features`
 *    correctly per side. `[datasetId]` is correctly `selectedDataset?.id` on the main instance,
 *    so this is an isolated copy-paste slip on `[features]` alone. Currently inert in practice:
 *    neither table declares `[actionColumns]`, and FilterableTableComponent only reads `features`
 *    for actionColumn-triggered fetchAndUpdate/geneSelected handling (confirmed in
 *    filterable-table.component.ts) — never to decide what the table itself renders, which reads
 *    only `data` — but it would silently make any future action column added to these tables
 *    write into the *compare* dataset's feature array while displaying the *main* dataset's
 *    stats. This component's `features` @Input() is a single generic per-instance value (the
 *    usual convention for this refactor), so whoever wires the main instance should bind
 *    `[features]="features"`, not `[features]="compareFeatures"`, to actually fix this — wiring
 *    this component in with the same literal expression the original template used would
 *    silently preserve the bug.
 * 2. FOUND, NOT FIXED (would require editing HexagonPlotComponent's fields/methods, out of scope):
 *    the TF Gene Network Graph and GRN Evaluation Plots panels rendered by this tab share their
 *    *state* with the "GRN Evaluation - On Demand" tab, not just their rendering code:
 *    - grnGraphFilter(Compare)/grnGraphDataPrior/Extended/Full(Compare) are written by BOTH this
 *      tab's loadPrecomputedGrnGraph (on tab open) AND the On-Demand tab's loadGrnGraph (on
 *      selecting a run) — selecting an on-demand run overwrites whatever this tab had cached, and
 *      the two tabs' filter-mode selectors both read/write the very same field.
 *    - grnPlots(Compare)/isLoadingGrnPlots(Compare) are likewise written by both
 *      loadPrecomputedGrnPlots (this tab) and loadGrnPlots (On-Demand tab).
 *    - Worse: visualizeGrnGraph always computes its target container id as
 *      `#tf_gene_net_container_<filterMode>(_compare)` — i.e. always one of *this tab's* three
 *      containers — even when invoked from the On-Demand tab's loadGrnGraph/
 *      visualizeGrnGraphForFilter call chain. The On-Demand tab's own graph container
 *      (`#grn_gene_net_container(_compare)`) is never targeted by any code, and that tab's
 *      `isLoadingGrnGraph(Compare)` flag is never set `true` anywhere either — see
 *      GrnEvaluationOnDemandPanelComponent's class doc for full detail. Practically: selecting an
 *      on-demand run's graph re-renders into *this* tab's TF Gene Network Graph containers
 *      (visible if this tab's DOM is still mounted, which `preserveContent="true"` on the
 *      enclosing mat-tab-group makes likely), not into the On-Demand tab the user is looking at.
 *    This component's grnGraphFilter/isLoadingGrnGraphPrior/Extended/Full/grnPlots/
 *    isLoadingGrnPlots @Input()s are independent, generic per-instance values (not hard-wired to
 *    assume sharing with a sibling On-Demand instance), so the wiring step CAN separate this state
 *    if desired, but wiring both tabs to today's same underlying fields will simply preserve the
 *    existing behavior.
 * 3. FOUND, NOT FIXED (informational only — symmetric on both sides, not a main/compare bug; would
 *    require editing HexagonPlotComponent): onChangeGrnGraphFilter's "no cached data for this
 *    filter mode yet" branch (the final `else` in its body) computes a `runName` local and never
 *    uses it — no fetch is triggered, so switching to an as-yet-unloaded filter mode silently
 *    shows stale/no data instead of loading it. This rarely bites *this* tab in practice
 *    (loadPrecomputedGrnGraph prefetches all 3 modes in parallel on tab-open, and its own
 *    `.then()` visualizes the current filter once its fetch resolves), but this component's
 *    `grnGraphFilterChange` output, once wired to the parent's onChangeGrnGraphFilter(value,
 *    isCompare), inherits this latent gap as-is.
 */
import { MatIconModule } from '@angular/material/icon';
import { MatButtonModule } from '@angular/material/button';
import { MatTooltipModule } from '@angular/material/tooltip';
import { InfoService } from '../info.service';

@Component({
  selector: 'app-grn-evaluation-panel',
  standalone: true,
  imports: [
    CommonModule,
    FilterableTableComponent,
    MatExpansionModule,
    MatFormField,
    MatLabel,
    MatOption,
    MatSelect,
    MatProgressSpinnerModule,
    MatIconModule,
    MatButtonModule,
    MatTooltipModule,
  ],
  templateUrl: './grn-evaluation-panel.component.html',
  styleUrls: ['./grn-evaluation-panel.component.scss'],
})
export class GrnEvaluationPanelComponent {
  constructor(public infoService: InfoService) {}

  @Input() isCompare = false;

  /** meta['peak_stats']?.[selectedDataset?.grn_evaluation_name || 'GRN_Evaluation'] (or
   * metaCompare/selectedDatasetCompare on the compare side) — the parent is expected to resolve
   * the eval-name lookup before passing this down, same as the original inline template
   * expression did. */
  @Input() peakStats: GrnTableData = [];
  /** !!meta['motif_stats'] (or metaCompare) — gates the Motif Statistics panel's *ngIf,
   * independent of whether motifStats itself resolves to data for the current eval name (mirrors
   * the original template's two separate expressions for panel visibility vs. table data). */
  @Input() hasMotifStats = false;
  @Input() motifStats: GrnTableData = [];

  /** See bug (1) above — the parent is expected to pass `features` on the main instance and
   * `compareFeatures` on the compare instance. */
  @Input() features: CellFeature[] = [];
  @Input() datasetId?: string;

  @Output() tableLoadingChange = new EventEmitter<{ loading: boolean; isCompare: boolean }>();
  /** Relayed from both `<app-table>`s' `featuresUpdated`. Parent binds
   * updateHexColors('#hexbin') (main) / updateHexColors('#hexbin-compare') (compare) at the two
   * call sites, matching the original template's pre-existing per-side handler. */
  @Output() featuresUpdated = new EventEmitter<void>();

  @Input() grnGraphFilter: GrnGraphFilterMode = 'full';
  /** Parent is expected to call its own onChangeGrnGraphFilter(value, isCompare) in response —
   * see bugs (2)/(3) above for what that currently does and doesn't do. */
  @Output() grnGraphFilterChange = new EventEmitter<GrnGraphFilterMode>();
  @Input() isLoadingGrnGraphPrior = false;
  @Input() isLoadingGrnGraphExtended = false;
  @Input() isLoadingGrnGraphFull = false;

  @Input() isLoadingGrnPlots = false;
  @Input() grnPlots: { [directory: string]: string[] } = {};
  /** Parent is expected to bind `(name) => this.getGrnImageUrl(name, isCompare)`. */
  @Input() getImageUrl: (imageName: string) => string = () => '';

  get priorContainerId(): string {
    return this.isCompare ? 'tf_gene_net_container_prior_compare' : 'tf_gene_net_container_prior';
  }

  get extendedContainerId(): string {
    return this.isCompare ? 'tf_gene_net_container_extended_compare' : 'tf_gene_net_container_extended';
  }

  get fullContainerId(): string {
    return this.isCompare ? 'tf_gene_net_container_full_compare' : 'tf_gene_net_container_full';
  }
}
