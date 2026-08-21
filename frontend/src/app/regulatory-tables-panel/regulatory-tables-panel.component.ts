import { Component, EventEmitter, Input, Output } from '@angular/core';
import { CommonModule } from '@angular/common';
import { FilterableTableComponent } from '../filterable-table/filterable-table.component';
import { CellFeature } from '../hexagon-view/cell-feature.types';

export type RegulatoryTablesPanelKind =
  | 'ligand_receptor'
  | 'cell_comp_tf_activity'
  | 'tf_activity'
  | 'pathway_activity';

type TableData = { [col: string]: { [index: string]: string | number } } | string[];

interface RegulatoryTableConfig {
  /** Key into the `meta`/`metaCompare` blob holding this table's rows. */
  metaKey: string;
  actionColumns: string[];
  emptyMessage: string;
  heading: string;
  /** Key into `selectedItemMap` for the "current selection" readout, and the `view` string
   * forwarded on `itemSelected` — independent of `metaKey` and `updateColumn` (e.g. the TF
   * Activity table reads from meta['tf_names'] but reports/selects under view 'tf_activity'). */
  view: string;
}

interface RegulatoryTablesPanelConfig {
  /** feature-property key FilterableTableComponent.fetchAndUpdate writes fetched scores under;
   * must match the selectedView/selectedCompareView string the color-by renderer reads. */
  updateColumn: string;
  primary: RegulatoryTableConfig;
  /** Only 'ligand_receptor' has a second table (NMF factors), toggled against the primary one. */
  secondary?: RegulatoryTableConfig;
  /** Static "Score: ..." text. Omitted for 'ligand_receptor', whose score line is dynamic
   * (driven by showGlobalScores) and built separately. */
  scoreLabel?: string;
  /** Whether the compare instance's heading gets " Comparison" appended — true for tf_activity/
   * pathway_activity, false for ligand_receptor/cell_comp_tf_activity, matching the original
   * markup's (inconsistent, but faithfully preserved) main/compare heading text. */
  appendComparisonSuffix: boolean;
}

const CONFIG: Record<RegulatoryTablesPanelKind, RegulatoryTablesPanelConfig> = {
  ligand_receptor: {
    updateColumn: 'ligand_receptor_relationships',
    primary: {
      metaKey: 'ligand_receptor_global_scores',
      actionColumns: ['ligand_receptor_cosine_similarity', 'ligand_receptor_p_value', 'ligand_receptor_category'],
      emptyMessage: 'No ligand-receptor global scores loaded',
      heading: 'Ligand–Receptor Global Scores',
      view: 'ligand_receptor_global_scores',
    },
    secondary: {
      metaKey: 'nmf_factors',
      actionColumns: ['ligand_receptor_NMF_factors'],
      emptyMessage: 'No ligand-receptor NMF factors loaded',
      heading: 'Ligand-Receptor NMF Scores',
      view: 'ligand_receptor_NMF_factors',
    },
    appendComparisonSuffix: false,
  },
  cell_comp_tf_activity: {
    updateColumn: 'cell_comp_tf_activity_similarity',
    primary: {
      metaKey: 'cell_comp_tf_activity_global_scores',
      actionColumns: ['cell_comp_tf_activity_cosine_similarity', 'cell_comp_tf_activity_category'],
      emptyMessage: 'No cell composition TF activity scores loaded',
      heading: 'Cell Composition TF Activity Global Scores',
      view: 'cell_comp_tf_activity_global_scores',
    },
    scoreLabel: 'Cell composition TF activity similarity',
    appendComparisonSuffix: false,
  },
  tf_activity: {
    updateColumn: 'tf_activity',
    primary: {
      metaKey: 'tf_names',
      actionColumns: ['tf_activity_score_ulm', 'tf_activity_padj_ulm'],
      emptyMessage: 'No TF activity scores loaded',
      heading: 'TF Activity',
      view: 'tf_activity',
    },
    scoreLabel: 'TF activity (ULM)',
    appendComparisonSuffix: true,
  },
  pathway_activity: {
    // NB: the original compare-side markup bound updateColumn to the typo 'pathwayActivity'
    // instead of 'pathway_activity' — see class doc below. Fixed here since this constant is
    // now shared by both instances.
    updateColumn: 'pathway_activity',
    primary: {
      metaKey: 'pathway_names',
      actionColumns: ['pathway_activity_score_mlm', 'pathway_activity_padj_mlm'],
      emptyMessage: 'No pathway activity scores loaded',
      heading: 'Pathway Activity',
      view: 'pathway_activity',
    },
    scoreLabel: 'Pathway activity (MLM)',
    appendComparisonSuffix: true,
  },
};

/**
 * Renders one of the four structurally-similar "regulatory scores" sidebar tabs — Ligand-Receptor
 * Relationships, Cell Composition TF Activity, TF Activity, and Pathway Activity — for one side
 * (main or compare). Extracted from HexagonPlotComponent, which used to hand-duplicate each tab's
 * markup verbatim for main vs. compare.
 *
 * A single `kind`-driven component was chosen over four separate ones: three of the four kinds
 * (cell_comp_tf_activity, tf_activity, pathway_activity) are byte-for-byte the same shape — one
 * "current selection" readout plus one `<app-table>` — differing only in which meta key/action
 * columns/updateColumn/copy they use. Only 'ligand_receptor' adds a second table (NMF factors)
 * behind a toggle. Centralizing the per-kind strings in one CONFIG map (rather than four templates
 * each repeating them) is also what surfaces the bugs below as a single fix instead of four.
 *
 * Bugs found in the original markup and fixed here (flagged per refactor convention, not silently):
 * 1. Functional bug: "Compare - Pathway Activity"'s `<app-table>` bound
 *    `[updateColumn]="'pathwayActivity'"` (camelCase) while every other instance — including the
 *    main-side Pathway Activity tab — used `'pathway_activity'`. FilterableTableComponent writes
 *    fetched scores to `feature.properties[updateColumn]`, and color-by rendering reads
 *    `feature.properties['pathway_activity']` (matching `selectedCompareView`), so selecting a
 *    pathway on the compare side wrote the score under a property the map never looks at,
 *    silently breaking compare-side pathway activity coloring. Fixed by using one shared
 *    `updateColumn` constant per kind.
 * 2. Cosmetic bug: every compare-side `<app-table>` across all four tabs omitted `[emptyMessage]`
 *    (present on every main-side instance), so an empty compare table showed a blank string
 *    instead of the explanatory message. Fixed by driving `emptyMessage` from the same per-kind
 *    config regardless of `isCompare`.
 *
 * Preserved as-is (not bugs, just pre-existing quirks kept faithfully):
 * - `showGlobalScores` (kind 'ligand_receptor' only) is a single boolean owned by
 *   HexagonPlotComponent, bound identically into both the main and compare instance — there was
 *   never a separate per-side toggle, so toggling either instance's button already flips both.
 *   Kept as an @Input()/@Output() pair rather than local state so that behavior is unchanged.
 * - Only tf_activity/pathway_activity append " Comparison" to the compare-side heading; the
 *   ligand_receptor and cell_comp_tf_activity headings are identical on both sides.
 * - The ligand_receptor "Selected Index"/"Selected index" readout label capitalization genuinely
 *   differs between the original main (capital I) and compare (lowercase i) markup; reproduced
 *   verbatim via `selectedIndexLabel` rather than normalized, since it has no functional effect.
 * - The ligand_receptor readout pipes each selected item through `json` before falling back with
 *   `||` (so a selection renders quoted, e.g. `"Slc2a1"`); the other three tabs don't. Preserved
 *   verbatim rather than unified, since changing it would change rendered text.
 */
@Component({
  selector: 'app-regulatory-tables-panel',
  standalone: true,
  imports: [CommonModule, FilterableTableComponent],
  templateUrl: './regulatory-tables-panel.component.html',
  styleUrls: ['./regulatory-tables-panel.component.scss'],
})
export class RegulatoryTablesPanelComponent {
  @Input({ required: true }) kind!: RegulatoryTablesPanelKind;
  @Input() isCompare = false;

  /** Raw `meta`/`metaCompare` blob; which key(s) are read from it is derived from `kind`. */
  @Input() meta: { [key: string]: any } = {};
  /** Raw `selectedItemByView`/`selectedItemByViewCompare` map, for the "current selection" readout. */
  @Input() selectedItemMap: { [view: string]: string | null } = {};
  @Input() features: CellFeature[] = [];
  @Input() datasetId?: string;

  @Input() showGlobalScores = true;
  @Output() showGlobalScoresChange = new EventEmitter<boolean>();

  /**
   * Emits in place of the original's per-call-site `onItemSelected($event, view[, true])`. The
   * parent is expected to call its own onItemSelected(event.event, event.view, isCompare) —
   * that also triggers fetchAndUpdate/color updates — not just record the selection itself.
   */
  @Output() itemSelected = new EventEmitter<{ event: { gene: string; action: string }; view: string }>();
  @Output() tableLoadingChange = new EventEmitter<{ loading: boolean; isCompare: boolean }>();
  /** Relayed from `<app-table>`; the parent binds updateHexColors() (main) / repaintBothViews()
   * (compare) at the two call sites, matching the original's pre-existing main/compare asymmetry. */
  @Output() featuresUpdated = new EventEmitter<void>();

  private get activeConfig(): RegulatoryTablesPanelConfig {
    return CONFIG[this.kind];
  }

  get hasSecondaryTable(): boolean {
    return !!this.activeConfig.secondary;
  }

  get selectionLabel(): string {
    return this.isCompare ? 'Comparison selection' : 'Current selection';
  }

  get selectedIndexLabel(): string {
    return this.kind === 'ligand_receptor' && !this.isCompare ? 'Selected Index' : 'Selected index';
  }

  get scoreLine(): string {
    return this.kind === 'ligand_receptor'
      ? (this.showGlobalScores ? 'Global scores' : 'NMF factors')
      : this.activeConfig.scoreLabel ?? '';
  }

  get updateColumn(): string {
    return this.activeConfig.updateColumn;
  }

  get primaryData(): TableData {
    return this.meta[this.activeConfig.primary.metaKey] ?? [];
  }

  get primaryActionColumns(): string[] {
    return this.activeConfig.primary.actionColumns;
  }

  get primaryEmptyMessage(): string {
    return this.activeConfig.primary.emptyMessage;
  }

  get primaryHeading(): string {
    return this.headingFor(this.activeConfig.primary);
  }

  get primarySelectedItem(): string | null {
    return this.selectedItemMap[this.activeConfig.primary.view];
  }

  get secondaryData(): TableData {
    const secondary = this.activeConfig.secondary;
    return (secondary ? this.meta[secondary.metaKey] : undefined) ?? [];
  }

  get secondaryActionColumns(): string[] {
    return this.activeConfig.secondary?.actionColumns ?? [];
  }

  get secondaryEmptyMessage(): string {
    return this.activeConfig.secondary?.emptyMessage ?? '';
  }

  get secondaryHeading(): string {
    const secondary = this.activeConfig.secondary;
    return secondary ? this.headingFor(secondary) : '';
  }

  get secondarySelectedItem(): string | null {
    const secondary = this.activeConfig.secondary;
    return secondary ? this.selectedItemMap[secondary.view] : null;
  }

  private headingFor(table: RegulatoryTableConfig): string {
    return this.isCompare && this.activeConfig.appendComparisonSuffix
      ? `${table.heading} Comparison`
      : table.heading;
  }

  toggleScoreView(): void {
    this.showGlobalScoresChange.emit(!this.showGlobalScores);
  }

  onPrimarySelected(event: { gene: string; action: string }): void {
    this.itemSelected.emit({ event, view: this.activeConfig.primary.view });
  }

  onSecondarySelected(event: { gene: string; action: string }): void {
    const secondary = this.activeConfig.secondary;
    if (secondary) this.itemSelected.emit({ event, view: secondary.view });
  }
}
