import { Component, EventEmitter, Input, Output } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';
import { SafeResourceUrl } from '@angular/platform-browser';
import { MatFormField, MatLabel } from '@angular/material/form-field';
import { MatOptgroup, MatOption } from '@angular/material/autocomplete';
import { MatSelect, MatSelectTrigger } from '@angular/material/select';
import { MatProgressSpinnerModule } from '@angular/material/progress-spinner';

/**
 * Renders the "Footprints" sidebar tab for one side (main or compare) — pre-computed footprint
 * plots, the on-demand "Compute New Footprint" form (cluster/group-by selector, a searchable
 * motif multi-select), and any resulting on-demand plots. Extracted from HexagonPlotComponent,
 * which used to hand-duplicate this markup verbatim for main vs. compare.
 *
 * The actual data-fetching stays in HexagonPlotComponent, not here:
 * - renderFootprintPlots (populates footprintPlotUrls/availableMotifs/availableCellTypes) is
 *   invoked from several cross-cutting call sites — the initial dataset-path subscription, cell
 *   clicks (displayCellDetails), cluster auto-selection (autoSelectDefaultCluster), and manual
 *   cluster selection (selectCluster) — matching the ClusterInfoPanelComponent/
 *   CoOccurrencePanelComponent precedent of keeping multi-call-site fetches in the parent.
 * - computeFootprint (the on-demand POST triggered by "Compute Footprints") also stays in the
 *   parent, for the same "actual fetch/compute of footprint plot URLs stays in the parent"
 *   reason. This component only emits `computeFootprintRequested` once the button is clicked;
 *   `footprintClusterBy`/`footprintMotifs` have already round-tripped to the parent via their
 *   `*Change` outputs by that point, so the parent's existing `computeFootprint(isCompare)` can
 *   run unmodified, reading its own fields exactly as before.
 *
 * `motifSearchQuery`/`footprintClusterBy`/`footprintMotifs` are @Input()/@Output() pairs — rather
 * than state owned locally by this component — specifically so the parent can keep resetting them
 * itself (as it already does on dataset reload, e.g. clearing footprintMotifs/motifSearchQuery)
 * simply by reassigning its own fields; the *Change outputs only relay user edits back up.
 *
 * KNOWN BUG, pre-existing in HexagonPlotComponent and NOT fixed here (out of this component's
 * scope — it only displays whatever @Input()s it is given):
 * Both the main and compare "Footprints" tabs currently bind their pre-computed plots to the
 * single `footprintPlotUrls` field. `footprintPlotUrlsCompare` exists on HexagonPlotComponent but
 * is never assigned — `renderFootprintPlots` always writes to `footprintPlotUrls`, regardless of
 * whether it was called with the main or the compare dataset — and `availableMotifs`/
 * `availableCellTypes`/`motifSearchQuery` are single fields shared by both sides rather than
 * per-side. In practice, whichever side last called renderFootprintPlots (e.g. a compare-side
 * cluster click) overwrites what the main tab displays too. This component takes single,
 * generically-named @Input()s per instance (same convention as CoOccurrencePanelComponent) so the
 * parent CAN supply correct per-side values once that source bug is fixed; wiring this component
 * in without also fixing renderFootprintPlots/splitting those fields will simply preserve the
 * existing cross-contamination.
 */
@Component({
  selector: 'app-footprint-panel',
  standalone: true,
  imports: [
    CommonModule,
    FormsModule,
    MatFormField,
    MatLabel,
    MatOptgroup,
    MatOption,
    MatSelect,
    MatSelectTrigger,
    MatProgressSpinnerModule,
  ],
  templateUrl: './footprint-panel.component.html',
  styleUrls: ['./footprint-panel.component.scss'],
})
export class FootprintPanelComponent {
  @Input() isCompare = false;

  /** Pre-computed footprint plot URLs (iframe srcs) for this dataset. */
  @Input() footprintPlotUrls: SafeResourceUrl[] = [];
  /** Plot URLs produced by this instance's on-demand "Compute Footprints" requests so far. */
  @Input() onDemandFootprintUrls: SafeResourceUrl[] = [];

  /** All motif ids available for this dataset, to populate the searchable multi-select. */
  @Input() availableMotifs: string[] = [];
  /** All cell types available for this dataset, to populate the "Specific Cell Type (is_…)" group. */
  @Input() availableCellTypes: string[] = [];

  @Input() motifSearchQuery = '';
  @Output() motifSearchQueryChange = new EventEmitter<string>();

  @Input() footprintClusterBy = 'cell_type';
  @Output() footprintClusterByChange = new EventEmitter<string>();

  @Input() footprintMotifs: string[] = [];
  @Output() footprintMotifsChange = new EventEmitter<string[]>();

  @Input() isComputingFootprint = false;
  @Input() footprintComputeError = '';

  /**
   * Emits when the user clicks "Compute Footprints". By this point footprintClusterBy/
   * footprintMotifs have already reached the parent via their `*Change` outputs, so it can call
   * its own computeFootprint(isCompare) directly in response — this component never performs the
   * HTTP request itself.
   */
  @Output() computeFootprintRequested = new EventEmitter<void>();

  get filteredMotifs(): string[] {
    const q = this.motifSearchQuery.trim().toLowerCase();
    return q ? this.availableMotifs.filter((m) => m.toLowerCase().includes(q)) : this.availableMotifs;
  }

  onMotifSearchInput(value: string): void {
    this.motifSearchQuery = value;
    this.motifSearchQueryChange.emit(value);
  }

  onClusterByChange(value: string): void {
    this.footprintClusterBy = value;
    this.footprintClusterByChange.emit(value);
  }

  onMotifsChange(value: string[]): void {
    this.footprintMotifs = value;
    this.footprintMotifsChange.emit(value);
  }
}
