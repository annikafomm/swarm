import { Component, Input } from '@angular/core';
import { CommonModule } from '@angular/common';
import { CellFeature, PropertyGroup, PropertyGroupItem } from '../hexagon-view/cell-feature.types';
import { InfoService } from '../info.service';
import { TranslationService } from '../translation.service';
import { MatIconModule } from '@angular/material/icon';
import { MatButtonModule } from '@angular/material/button';
import { MatTooltipModule } from '@angular/material/tooltip';
import { MatExpansionModule } from '@angular/material/expansion';

/**
 * Renders the "Cell Information" sidebar tab for one side (main or compare). Extracted from
 * HexagonPlotComponent, which used to share this via a single ng-template + ngTemplateOutlet +
 * context object (one class, two call sites) — now each instance just takes its own `cell`/
 * `currentView` directly. Lives in the sidebar (`.metadata-tables-container`), a separate DOM
 * region from the map (`#comparison-box`) — NOT nested inside HexagonViewComponent, which only
 * covers the map itself. See the refactor plan's "CORRECTION" section for why.
 */
@Component({
  selector: 'app-cell-info-panel',
  standalone: true,
  imports: [CommonModule, MatIconModule, MatButtonModule, MatTooltipModule, MatExpansionModule],
  templateUrl: './cell-info-panel.component.html',
  styleUrls: ['./cell-info-panel.component.scss'],
})
export class CellInfoPanelComponent {
  @Input() cell: CellFeature | null = null;
  @Input() currentView = '';
  @Input() datasetFeatures: CellFeature[] = [];

  constructor(
    public infoService: InfoService,
    private translationService: TranslationService,
  ) { }

  private readonly hiddenPropKeys = new Set<string>(['barcode', 'cell_type', 'leiden', 'centroid', 'color']);
  private expandedProps = new Set<string>();
  private expandedGroupKeys = new Set<string>();
  private lastCellBarcode: string | null = null;

  // Cache non-empty properties per datasetFeatures reference
  private nonEmptyPropsCache = new Set<string>();
  private lastFeaturesRef: CellFeature[] | null = null;

  shouldShowProperty(key: string): boolean {
    if (key == null) return true;
    const k = String(key).toLowerCase();
    return !this.hiddenPropKeys.has(k);
  }

  public isValueEmpty(v: unknown): boolean {
    if (v === null || v === undefined) return true;
    if (typeof v === 'string') {
      const trimmed = v.trim();
      return (
        trimmed === '' ||
        trimmed.toLowerCase() === 'nan' ||
        trimmed.toLowerCase() === 'null' ||
        trimmed.toLowerCase() === 'none' ||
        trimmed.toLowerCase() === 'undefined'
      );
    }
    if (typeof v === 'number') {
      return !Number.isFinite(v);
    }
    if (Array.isArray(v)) {
      return v.length === 0 || v.every((x) => this.isValueEmpty(x));
    }
    if (typeof v === 'object') {
      return Object.keys(v).length === 0;
    }
    return false;
  }

  private computeNonEmptyProps(): void {
    if (this.datasetFeatures === this.lastFeaturesRef) return;
    this.lastFeaturesRef = this.datasetFeatures;
    this.nonEmptyPropsCache.clear();

    if (!this.datasetFeatures || this.datasetFeatures.length === 0) return;

    // Collect all property keys that appear anywhere in the dataset
    const allKeys = new Set<string>();
    for (let i = 0; i < this.datasetFeatures.length; i++) {
      const p = this.datasetFeatures[i]?.properties;
      if (p) {
        for (const k of Object.keys(p)) {
          allKeys.add(k);
        }
      }
    }

    // A property is considered present in the dataset if at least one cell has a non-empty value
    for (const key of allKeys) {
      const hasValue = this.datasetFeatures.some((f) => {
        const val = f.properties?.[key];
        return !this.isValueEmpty(val);
      });
      if (hasValue) {
        this.nonEmptyPropsCache.add(key);
      }
    }
  }

  /**
   * Buckets a cell's raw properties into biologist-facing sections. Not memoized: fetchAndUpdate
   * (still in HexagonPlotComponent) mutates `cell.properties` in place on the already-selected
   * cell (e.g. when switching the colored-by gene/score without deselecting the cell), so caching
   * by cell identity would show stale values for the "Current View Value" section.
   */
  groupedCellProperties(cell: CellFeature | null): PropertyGroup[] {
    if (!cell) return [];
    const groups = this.computeGroupedProperties(cell);
    const barcode = cell.properties.barcode;
    if (barcode !== this.lastCellBarcode) {
      this.lastCellBarcode = barcode;
      this.expandedGroupKeys.clear();
      // Default-expand first 3 groups for newly selected cell
      for (let i = 0; i < Math.min(3, groups.length); i++) {
        this.expandedGroupKeys.add(groups[i].key);
      }
    }
    return groups;
  }

  isGroupExpanded(key: string): boolean {
    return this.expandedGroupKeys.has(key);
  }

  onGroupOpened(key: string): void {
    this.expandedGroupKeys.add(key);
  }

  onGroupClosed(key: string): void {
    this.expandedGroupKeys.delete(key);
  }

  trackByGroup(index: number, group: PropertyGroup): string {
    return group.key;
  }

  trackByItem(index: number, item: PropertyGroupItem): string {
    return item.key;
  }

  private computeGroupedProperties(cell: CellFeature): PropertyGroup[] {
    this.computeNonEmptyProps();

    const buckets = new Map<string, PropertyGroupItem[]>();
    for (const key of Object.keys(cell.properties)) {
      if (!this.shouldShowProperty(key)) continue;

      // If datasetFeatures is provided and this property is 100% empty across the entire dataset, skip it
      if (this.datasetFeatures.length > 0 && !this.nonEmptyPropsCache.has(key)) {
        continue;
      }

      const rawValue = (cell.properties as any)[key];
      const category = this.infoService.getCategoryForKey(key, rawValue);
      const item: PropertyGroupItem = {
        key,
        label: this.label(key),
        value: rawValue,
        info: this.infoService.getFieldInfo(key),
      };
      const bucket = buckets.get(category);
      if (bucket) {
        bucket.push(item);
      } else {
        buckets.set(category, [item]);
      }
    }

    return this.infoService.categories
      .map((def) => ({
        ...def,
        items: (buckets.get(def.key) ?? []).sort((a, b) => a.label.localeCompare(b.label, 'de', { sensitivity: 'base' })),
      }))
      .filter((group) => group.items.length > 0);
  }

  // Also used by HexagonPlotComponent's renderLegend, which has its own copy — kept as a tiny
  // duplicated pure function rather than shared, to avoid a cross-file dependency for six lines.
  label(key: string): string {
    const translated = this.translationService.translateSync(key);
    if (translated !== key) return translated;
    return key.replace(/_/g, ' ').replace(/\b\w/g, c => c.toUpperCase());
  }

  // ----- Property-value formatting helpers (used only by the propertyValue ng-template) -----

  isArray(v: any): v is any[] {
    return Array.isArray(v);
  }

  isNestedArray(v: any): v is any[][] {
    return Array.isArray(v) && v.length > 0 && v.every(row => Array.isArray(row) || this.looksLikeArrayString(row));
  }

  isNumberLike(v: unknown): v is number | string {
    return (typeof v === 'number' && Number.isFinite(v)) ||
      (typeof v === 'string' && v.trim() !== '' && Number.isFinite(+v));
  }

  toNumberLike(v: number | string): number {
    return typeof v === 'number' ? v : Number(v);
  }

  isNumericArray(arr: any): arr is (number | string)[] {
    return Array.isArray(arr) && arr.length > 0 && arr.every(x => this.isNumberLike(x));
  }

  isPrimitive(v: unknown): v is string | number | boolean | null {
    return v === null || ['string', 'number', 'boolean'].includes(typeof v as string);
  }

  getArrayStats(arr: (number | string)[]) {
    const nums = arr.map(x => Number(this.toNumberLike(x)));
    const min = Math.min(...nums);
    const max = Math.max(...nums);
    const avg = nums.reduce((a, b) => a + b, 0) / nums.length;
    return { min, max, avg };
  }

  formatValue(v: unknown): string {
    if (Array.isArray(v)) return v.join(', ');
    if (v && typeof v === 'object') return JSON.stringify(v, null, 2);
    return String(v);
  }

  toggleExpand(key: string) {
    if (this.expandedProps.has(key)) this.expandedProps.delete(key);
    else this.expandedProps.add(key);
  }

  isExpanded(key: string) {
    return this.expandedProps.has(key);
  }

  trackByIndex(index: number) { return index; }

  asArrayRow(row: any): any[] {
    if (Array.isArray(row)) return row;
    if (this.looksLikeArrayString(row)) {
      try {
        const parsed = JSON.parse(row as string);
        return Array.isArray(parsed) ? parsed : [row];
      } catch {
        return [row];
      }
    }
    return [row];
  }

  private looksLikeArrayString(v: any): v is string {
    return typeof v === 'string' && /^\s*\[.*\]\s*$/.test(v);
  }

  toJsonCompact(obj: unknown, max = 120): string {
    try {
      const s = JSON.stringify(obj);
      return s.length > max ? s.slice(0, max) + '…' : s;
    } catch {
      return String(obj);
    }
  }

  isPlainObject(v: any): v is Record<string, any> {
    return v !== null && typeof v === 'object' && !Array.isArray(v);
  }

  objectKeyCount(obj: Record<string, any>): number {
    return Object.keys(obj).length;
  }

  objectEntries(obj: Record<string, any>): Array<{ key: string; value: any }> {
    return Object.keys(obj).sort().map(k => ({ key: k, value: obj[k] }));
  }

  prettyKey(k: string): string {
    return k.replace(/[_-]+/g, ' ')
      .replace(/\b\w/g, c => c.toUpperCase());
  }

  dictId(propLabel: string): string {
    return `DICT::${propLabel}`;
  }
}
