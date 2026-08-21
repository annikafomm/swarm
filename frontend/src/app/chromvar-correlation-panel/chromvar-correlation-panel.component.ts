import { Component, EventEmitter, Input, Output } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FilterableTableComponent } from '../filterable-table/filterable-table.component';
import { CellFeature } from '../hexagon-view/cell-feature.types';

type TableData = { [col: string]: { [index: string]: string | number } } | string[];

/**
 * Renders the "ChromVar spatial correlation : Moran's I / Geary's C" sidebar tab for one side
 * (main or compare) — a toggle button switching between the two tables, each scored under the
 * same `chromvar_total_sum` color-by property. Extracted from HexagonPlotComponent, which used
 * to hand-duplicate this markup verbatim for main vs. compare.
 *
 * The Moran's I / Geary's C toggle (`showChromvarMoranI`/`compareShowChromvarMoranI` in the
 * original) is moved to purely local component state rather than kept as an @Input()/@Output()
 * pair: unlike RegulatoryTablesPanelComponent's `showGlobalScores` (one boolean shared identically
 * by both the main and compare instance), this tab's main and compare booleans were already two
 * independent fields, read nowhere outside this tab's own markup — genuinely self-contained
 * per-instance UI state with nothing left for the parent to orchestrate.
 */
import { MatIconModule } from '@angular/material/icon';
import { MatButtonModule } from '@angular/material/button';
import { MatTooltipModule } from '@angular/material/tooltip';
import { InfoService } from '../info.service';

@Component({
  selector: 'app-chromvar-correlation-panel',
  standalone: true,
  imports: [CommonModule, FilterableTableComponent, MatIconModule, MatButtonModule, MatTooltipModule],
  templateUrl: './chromvar-correlation-panel.component.html',
  styleUrls: ['./chromvar-correlation-panel.component.scss'],
})
export class ChromvarCorrelationPanelComponent {
  constructor(public infoService: InfoService) {}

  @Input() isCompare = false;
  /** `meta['chromvar_moranI']` / `metaCompare['chromvar_moranI']`. */
  @Input() moranIData: TableData = [];
  /** `meta['chromvar_gearyC']` / `metaCompare['chromvar_gearyC']`. */
  @Input() gearyCData: TableData = [];
  @Input() features: CellFeature[] = [];
  @Input() datasetId?: string;

  /**
   * Relayed from `<app-table>`; the parent is expected to bind updateHexColors() (main) /
   * updateHexColors('#hexbin-compare') (compare) at the two call sites, matching the original's
   * pre-existing main/compare asymmetry.
   */
  @Output() featuresUpdated = new EventEmitter<void>();
  @Output() tableLoadingChange = new EventEmitter<{ loading: boolean; isCompare: boolean }>();

  showMoranI = true;

  toggleView(): void {
    this.showMoranI = !this.showMoranI;
  }
}
