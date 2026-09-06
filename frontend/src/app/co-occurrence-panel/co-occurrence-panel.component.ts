import { Component, EventEmitter, Input, Output } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FormsModule } from '@angular/forms';

/**
 * Renders the co-occurrence matrix embedded in the "Cluster Information" sidebar tab, for one
 * side (main or compare) — interval selector, summary stats, and the per-cluster co-occurrence
 * table. Extracted from HexagonPlotComponent, which used to hand-duplicate this markup verbatim
 * for main vs. compare.
 *
 * Like ClusterInfoPanelComponent, this is a genuinely passive, @Input()-driven display — fetching
 * the co-occurrence matrix itself (updateCoOccurrenceTable) stays in HexagonPlotComponent, because
 * it's also invoked from cell clicks and cluster selection, both of which touch cross-cutting
 * state (GRN subgraphs, footprint plots, map highlighting) that belongs to later refactor phases.
 * Only the pure computation that reads solely from this component's own @Input()s
 * (getCoOccurrenceColor, getIntervalStats, getActualIntervalRange/formatActualIntervalRange) moves
 * here, alongside the interval navigation, which only affects which interval this tab displays.
 *
 * Bug fix vs. the original: the compare-side block used the *main* side's `clusterCount` for its
 * table header (`<th>` loop), and its interval dropdown/navigation used an off-by-one range that
 * disagreed with the compare-side's own (correct) one — see `maxInterval`/`clusterCount` below.
 * Taking a single value per @Input() per instance means the parent must supply the correct
 * per-side value, and both sides now render with the same (corrected) navigation logic.
 */
import { MatIconModule } from '@angular/material/icon';
import { MatButtonModule } from '@angular/material/button';
import { MatTooltipModule } from '@angular/material/tooltip';
import { InfoService } from '../info.service';

@Component({
  selector: 'app-co-occurrence-panel',
  standalone: true,
  imports: [CommonModule, FormsModule, MatIconModule, MatButtonModule, MatTooltipModule],
  templateUrl: './co-occurrence-panel.component.html',
  styleUrls: ['./co-occurrence-panel.component.scss'],
})
export class CoOccurrencePanelComponent {
  constructor(public infoService: InfoService) {}

  @Input() isCompare = false;
  /** The cluster whose co-occurrence row this instance displays. */
  @Input() selectedClusterId: string | number | null = null;
  /** Co-occurrence value of `selectedClusterId` against each target cluster, at `selectedInterval`. */
  @Input() coOccurrenceData: number[] = [];
  /** Number of target-cluster columns to render in the table header. */
  @Input() clusterCount = 10;
  /** 75th-percentile value (computed by the parent) above which a cell is highlighted. */
  @Input() threshold = 0.5;
  /** Highest valid interval index for this dataset. */
  @Input() maxInterval = 49;
  @Input() selectedInterval = 0;
  /** Cleaned (numeric, finite) interval boundary values for this dataset, if any. */
  @Input() intervalBoundaries: number[] = [];
  /** Continuous color scale (already domained to [0, 1]) used to shade table cells by intensity. */
  @Input() continuousColorScale!: (value: number) => string;

  /**
   * Emits the new interval index when the user changes the dropdown or clicks Previous/Next. The
   * parent is expected to assign it to its own selectedInterval/selectedIntervalCompare and call
   * updateCoOccurrenceTable(isCompare) in response — this component stays passive and does not
   * refetch data itself.
   */
  @Output() intervalChanged = new EventEmitter<number>();

  /**
   * Emits when the user clicks a table cell to jump to that target cluster. The parent is expected
   * to call its own selectCluster(id, isCompare) — same convention as ClusterInfoPanelComponent's
   * clusterSelected.
   */
  @Output() clusterSelected = new EventEmitter<number>();

  get intervalOptions(): number[] {
    return Array.from({ length: this.maxInterval + 1 }, (_, i) => i);
  }

  get clusterIndices(): number[] {
    return Array.from({ length: this.clusterCount }, (_, i) => i);
  }

  selectPrevInterval(): void {
    if (this.selectedInterval > 0) {
      this.intervalChanged.emit(this.selectedInterval - 1);
    }
  }

  selectNextInterval(): void {
    if (this.selectedInterval < this.maxInterval) {
      this.intervalChanged.emit(this.selectedInterval + 1);
    }
  }

  getIntervalStats(): { min: number; max: number; avg: number } {
    const allValues = this.coOccurrenceData.flat().filter((val) => val > 0);
    if (allValues.length === 0) return { min: 0, max: 0, avg: 0 };

    const min = Math.min(...allValues);
    const max = Math.max(...allValues);
    const avg = allValues.reduce((sum, val) => sum + val, 0) / allValues.length;

    return { min, max, avg: Math.round(avg * 100) / 100 };
  }

  getCoOccurrenceColor(value: number): string {
    if (value === 0) return '#f8f9fa';

    const coOccurrenceWithoutSameCluster = this.coOccurrenceData.filter(
      (_, index) => index !== this.selectedClusterId,
    );
    const maxValue = Math.max(...coOccurrenceWithoutSameCluster);
    const intensity = Math.min(value / maxValue, 1);
    return this.continuousColorScale(intensity);
  }

  getActualIntervalRange(index: number): { start: number; end: number } | null {
    const intervals = this.intervalBoundaries;
    if (!intervals.length || index < 0 || index >= intervals.length) return null;

    return {
      start: index === 0 ? 0 : intervals[index - 1],
      end: intervals[index],
    };
  }

  formatActualIntervalRange(index: number, digits = 1): string {
    const range = this.getActualIntervalRange(index);
    if (!range) return '';

    return `(${range.start.toFixed(digits)}–${range.end.toFixed(digits)})`;
  }

  tooltipFor(targetCluster: number, value: number): string {
    return `Cluster ${this.selectedClusterId} → Cluster ${targetCluster}: ${value} (Click to select)`;
  }
}
