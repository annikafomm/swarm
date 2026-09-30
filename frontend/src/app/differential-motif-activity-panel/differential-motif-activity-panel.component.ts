import { Component, Input } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FilterableTableComponent } from '../filterable-table/filterable-table.component';

type TableData = { [col: string]: { [index: string]: string | number } } | string[];

/**
 * Renders the "Differential Motif Activity" sidebar tab for one side (main or compare) — one
 * heading + read-only `<app-table>` per comparison group. Extracted from HexagonPlotComponent,
 * which used to hand-duplicate this markup verbatim for main vs. compare.
 *
 * Kept as its own component rather than folded into ChromvarCorrelationPanelComponent: despite
 * both tabs living under the ChromVar/motif umbrella, their shapes share nothing — this one has
 * no toggle state and iterates a dynamic, unbounded set of comparisons (`| keyvalue`) instead of
 * two fixed named tables, and its `<app-table>` doesn't bind `features`/`updateColumn`/
 * `datasetId`/`isCompare` or any output at all (no `actionColumns`, so nothing on these tables is
 * ever actionable) — identically on both main and compare, so this is a pre-existing read-only
 * display, not a copy-paste omission. Preserved verbatim.
 *
 * The original's `asTableData()` cast helper on HexagonPlotComponent (needed only because
 * `meta`/`metaCompare` are untyped `{ [key: string]: any }` blobs) is unnecessary here since
 * `topMotifs` is typed directly.
 */
@Component({
  selector: 'app-differential-motif-activity-panel',
  standalone: true,
  imports: [CommonModule, FilterableTableComponent],
  templateUrl: './differential-motif-activity-panel.component.html',
  styleUrls: ['./differential-motif-activity-panel.component.scss'],
})
export class DifferentialMotifActivityPanelComponent {
  /** `meta['diff_motif_activity_top_motifs']` / `metaCompare['diff_motif_activity_top_motifs']`. */
  @Input() topMotifs: { [comparisonKey: string]: TableData } | null | undefined = null;
}
