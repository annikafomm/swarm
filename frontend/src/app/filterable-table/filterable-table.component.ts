import { Component, OnInit, Input, Output, EventEmitter, OnChanges, SimpleChanges, ElementRef } from '@angular/core';
import { FormsModule } from '@angular/forms';
import { CommonModule } from '@angular/common';
import { SessionService } from '../session.service';
import { HttpClient } from '@angular/common/http';
import { TranslatePipe } from '../translate.pipe';

import { MatIconModule } from '@angular/material/icon';
import { MatButtonModule } from '@angular/material/button';
import { MatTooltipModule } from '@angular/material/tooltip';
import { InfoService } from '../info.service';

@Component({
  selector: 'app-table',
  templateUrl: './filterable-table.component.html',
  styleUrls: ['./filterable-table.component.scss'],
  standalone: true,
  imports: [FormsModule, CommonModule, TranslatePipe, MatIconModule, MatButtonModule, MatTooltipModule],
})
export class FilterableTableComponent implements OnInit, OnChanges {
  /**
   * If data is an object, treat as {col: {index: value}}
   * If data is a string array, treat as index column
   */
  @Input() data!:
    | { [col: string]: { [index: string]: string | number } }
    | string[];
  @Input() actionColumns: string[] = [];
  @Input() features!: any;
  @Input() updateColumn!: string;
  @Input() datasetId?: string;
  @Input() isLoading: boolean = false;
  @Input() emptyMessage: string = '';
  @Input() isCompare: boolean = false;
  /** Hides the auto-generated "Index" column/filter — for tables where the row index is a
   * meaningless post-merge row number rather than a searchable identifier (e.g. the ChromVAR
   * correlation tables, which already expose TF/motif_id as their own filterable columns). */
  @Input() showIndexColumn: boolean = true;
  /** Header text for the index column — defaults to the generic "Index" but should be set to
   * whatever the row index actually identifies for a given table (e.g. "Gene", "TF"). */
  @Input() indexColumnLabel: string = 'Index';
  /**
   * Filter/sort state, as plain objects rather than primitives specifically so two tables that
   * want to stay in sync (e.g. the ChromVAR Moran's I / Geary's C toggle, which swaps which
   * table is in the DOM via *ngIf rather than keeping both alive) can be handed the *same*
   * object reference from their parent — mutations from either table are then visible to both,
   * with no @Output()/two-way-binding plumbing needed. Tables that don't bind these get their
   * own private, independent object per the field initializer below (unchanged default behavior).
   */
  @Input() filters: { [col: string]: string } = {};
  @Input() sortState: { column: string | null; asc: boolean } = { column: null, asc: true };
  /** False for tables sharing state across a data swap that's really "the same rows, viewed
   * differently" (see above) — everyone else keeps the default (reset filters/sort whenever a
   * genuinely new dataset arrives). */
  @Input() resetStateOnDataChange: boolean = true;
  @Output() featuresUpdated = new EventEmitter<void>();
  @Output() geneSelected = new EventEmitter<{ gene: string; action: string }>();
  @Output() geneSelectedCompare = new EventEmitter<{ gene: string; action: string }>();
  @Output() loadingChange = new EventEmitter<{ loading: boolean; isCompare: boolean }>();
  constructor(
    private http: HttpClient,
    private sessionService: SessionService,
    private el: ElementRef,
    public infoService: InfoService,
  ) { }

  columns: string[] = [];
  rows: any[] = [];
  availableActionColumns: string[] = [];
  private readonly virtualActionColumns = new Set<string>([
    'gene_expression',
    'chromvar_spot_scores',
    'tf_activity_score_ulm',
    'tf_activity_padj_ulm',
    'pathway_activity_score_mlm',
    'pathway_activity_padj_mlm',
    'cell_comp_tf_activity_cosine_similarity',
    'cell_comp_tf_activity_category',
    'ligand_receptor_cosine_similarity',
    'ligand_receptor_p_value',
    'ligand_receptor_category',
    'ligand_receptor_NMF_factors',
    'ligand_receptor_cosine_similarity',
    'ligand_receptor_p_value',
    'ligand_receptor_category'
  ]);

  // Pagination
  pageSize = 6; // number of rows per page
  currentPage = 1;

  ngOnInit() {
    this.prepareTable();
  }

  ngOnChanges(changes: SimpleChanges) {
    // Rebuild table whenever the provided data or features change
    if (changes['data'] || changes['features']) {
      this.currentPage = 1;
      if (this.resetStateOnDataChange) {
        this.sortState.column = null;
        this.sortState.asc = true;
        this.filters = {};
      }
      this.prepareTable();
    }
  }

  prepareTable() {
    if (!this.data) {
      this.columns = [];
      this.rows = [];
      return;
    }
    else if (Array.isArray(this.data)) {
      this.columns = [];
      this.rows = this.data.map((val) => ({
        index: val,
      }));
    } else {
      const tableData = this.data as {
        [col: string]: { [index: string]: string | number };
      };

      this.columns = Object.keys(tableData);
      const indexes = Object.keys(tableData[this.columns[0]] || {});

      this.rows = indexes.map((idx) => {
        const row: any = { index: idx };
        this.columns.forEach((col) => {
          row[col] = tableData[col][idx];
        });
        return row;
      });
    }

    // Filter action columns to only show those that are available in the data
    this.filterAvailableActionColumns();
  }

  filterAvailableActionColumns() {
    if (!this.data || (Array.isArray(this.data) && this.data.length === 0)) {
      this.availableActionColumns = [];
      return;
    }

    if ('motif_id' in this.data) {
      this.availableActionColumns = this.actionColumns;
    } else if (!Array.isArray(this.data)) {
      const tableData = this.data as {
        [col: string]: { [index: string]: string | number };
      };
      const matchingColumns = this.actionColumns.filter(
        (col) => col in tableData
      );
      const virtualColumns = this.actionColumns.filter(
        (col) => this.virtualActionColumns.has(col)
      );

      // Show only actually available columns, plus known virtual actions.
      this.availableActionColumns = Array.from(
        new Set([...matchingColumns, ...virtualColumns])
      );
    } else {
      this.availableActionColumns = this.actionColumns;
    }


    // Check if action columns exist in the actual data
    // if (!Array.isArray(this.data)) {
    //   const tableData = this.data as {
    //     [col: string]: { [index: string]: string | number };
    //   };
    //   this.availableActionColumns = this.actionColumns.filter(
    //     (col) => col in tableData
    //   );
    // } else {
    //   this.availableActionColumns = this.actionColumns;
    // }


    if (Array.isArray(this.data)) {
      this.availableActionColumns = this.actionColumns;
      return;
    }

    const tableData = this.data as {
      [col: string]: { [index: string]: string | number };
    };

    this.availableActionColumns = this.actionColumns.filter((col) => {
      if (col === 'show_on_plot') return true;
      if (col === 'gene_expression') return true;
      if (col === 'chromvar_spot_scores') return true;
      if (this.virtualActionColumns.has(col)) return true;  // ADD THIS LINE

      return col in tableData;
    });
  }

  /**
   * Whether `col` should get a text-filter input. Checks every row, not just the first — a
   * column that's a string everywhere except a null/undefined value in row 0 (e.g. an
   * unannotated first entry) would otherwise silently lose its search box for all rows.
   */
  isSearchableColumn(col: string): boolean {
    return this.rows.some((row) => typeof row[col] === 'string');
  }

  hasData(): boolean {
    return this.rows && this.rows.length > 0;
  }

  displayColumnName(col: string): string {
    return col.replace(/_/g, ' ');
  }

  truncateMiddle(text: string, front = 4, back = 7): string {
    // Render long ensembl IDs as ENSG...78531.1
    if (!text.startsWith('ENS')) return text;
    return text.slice(0, front) + '…' + text.slice(text.length - back);
  }

  displayNumeric(value: any): string {
    if (typeof value === 'number') {
      if (Number.isNaN(value)) {
        return 'N/A';
      }
      if (value === 0) {
        return '0.0';
      }

      // Use scientific notation for very small or very large numbers
      if (Math.abs(value) < 0.001 || Math.abs(value) >= 1e5) {
        const exp = value.toExponential(2); // e.g., "2.34e-5"
        const [mantissa, exponent] = exp.split('e');
        const expNum = parseInt(exponent, 10);

        // Map digits and minus sign to Unicode superscripts
        const superscripts: Record<string, string> = {
          '-': '⁻',
          '0': '⁰',
          '1': '¹',
          '2': '²',
          '3': '³',
          '4': '⁴',
          '5': '⁵',
          '6': '⁶',
          '7': '⁷',
          '8': '⁸',
          '9': '⁹',
        };

        const expStr = expNum
          .toString()
          .split('')
          .map((ch) => superscripts[ch] ?? ch)
          .join('');

        return `${mantissa} × 10${expStr}`;
      }

      // Otherwise, fixed 3 decimals
      return value.toFixed(3);
    }
    return value;
  }

  get filteredRows() {
    let result = [...this.rows];

    // filtering (strings only)
    for (const col of ['index', ...this.columns]) {
      if (this.filters[col]) {
        result = result.filter(
          (row) =>
            typeof row[col] === 'string' &&
            row[col].toLowerCase().includes(this.filters[col].toLowerCase()),
        );
      }
    }

    // sorting
    if (this.sortState.column) {
      const sortKey = this.sortState.column;
      const asc = this.sortState.asc;
      result.sort((a, b) => {
        const valA = a[sortKey];
        const valB = b[sortKey];

        if (typeof valA === 'number' && typeof valB === 'number') {
          return asc ? valA - valB : valB - valA;
        }
        return asc
          ? String(valA).localeCompare(String(valB))
          : String(valB).localeCompare(String(valA));
      });
    }

    return result;
  }

  // Return only rows for current page
  get pagedRows() {
    const start = (this.currentPage - 1) * this.pageSize;
    return this.filteredRows.slice(start, start + this.pageSize);
  }

  get totalPages() {
    return Math.ceil(this.filteredRows.length / this.pageSize) || 1;
  }

  setPage(page: number) {
    if (page >= 1 && page <= this.totalPages) {
      this.currentPage = page;
    }
  }

  toggleSort(col: string) {
    if (this.sortState.column === col) {
      this.sortState.asc = !this.sortState.asc;
    } else {
      this.sortState.column = col;
      this.sortState.asc = true;
    }
  }

  setSort(col: string, ascending: boolean) {
    this.sortState.column = col;
    this.sortState.asc = ascending;
  }

  // async fetchAndUpdate(columnName: string, index: string) {
  //   let isGeneExpression = columnName === 'gene_expression';
  //   let request = isGeneExpression
  //     ? `${this.sessionService.apiUrl}/X/${index}`
  //     : `${this.sessionService.apiUrl}/obsm/${columnName}/${index}`;
  //   this.sessionService
  //     .callWithSession(() => this.http.get(request, { withCredentials: true }))
  //     .subscribe({
  //       next: (res) => {
  //         const data = res as { [barcode: string]: any };

  //         if (this.features) {
  //           for (const feature of this.features) {
  //             const barcode = feature.properties?.barcode;
  //             if (barcode && data[barcode] !== undefined) {
  //               feature.properties[this.updateColumn] = data[barcode];
  //             }
  //           }
  //         }
  //         if (isGeneExpression) {
  //           console.log(`[Backend] Loaded adata[:, ${index}].X`);
  //         } else {
  //           console.log(`[Backend] Loaded adata.obsm[${columnName}][${index}]`);
  //         }

  //         this.featuresUpdated.emit();
  //       },
  //       error: (err) => {
  //         if (isGeneExpression) {
  //           console.error(`[Backend] Failed to load adata[:, ${index}].X`, err);
  //         } else {
  //           console.error(
  //             `[Backend] Failed to load adata.obsm[${columnName}][${index}]`,
  //             err,
  //           );
  //         }
  //       },
  //     });
  // }

  async fetchAndUpdate(columnName: string, index: string) {
    const isGeneExpression = columnName === 'gene_expression';
    const isChromvar = columnName === 'chromvar_spot_scores';

    const safeIndex = encodeURIComponent(index);
    const datasetQuery = this.datasetId ? `?dataset_id=${encodeURIComponent(this.datasetId)}` : '';

    const request = isGeneExpression
      ? `${this.sessionService.apiUrl}/X/${safeIndex}${datasetQuery}`
      : isChromvar
        ? `${this.sessionService.apiUrl}/obsm/chromvar_spot_scores/${safeIndex}${datasetQuery}`
        : `${this.sessionService.apiUrl}/obsm/${encodeURIComponent(columnName)}/${safeIndex}${datasetQuery}`;

    this.loadingChange.emit({ loading: true, isCompare: this.isCompare });

    this.sessionService
      .callWithSession(() => this.http.get(request, { withCredentials: true }))
      .subscribe({
        next: (res) => {
          const data = res as { [barcode: string]: any };

          if (this.features) {
            for (const feature of this.features) {
              const barcode = feature.properties?.barcode;
              if (barcode && data[barcode] !== undefined) {
                feature.properties[this.updateColumn] = data[barcode];
              }
            }
          }

          if (isGeneExpression) {
            console.log(`[Backend] Loaded adata[:, ${index}].X`);
          } else if (isChromvar) {
            console.log(`[Backend] Loaded ChromVAR score for motif '${index}'`);
          } else {
            console.log(`[Backend] Loaded adata.obsm[${columnName}][${index}]`);
          }

          this.loadingChange.emit({ loading: false, isCompare: this.isCompare });
          this.featuresUpdated.emit();
        },
        error: (err) => {
          if (this.features) {
            for (const feature of this.features) {
              feature.properties[this.updateColumn] = undefined;
            }
            this.featuresUpdated.emit();
          }

          this.loadingChange.emit({ loading: false, isCompare: this.isCompare });

          if (isGeneExpression) {
            console.error(`[Backend] Failed to load adata[:, ${index}].X`, err);
          } else if (isChromvar) {
            console.error(`[Backend] Failed to load ChromVAR score for motif '${index}'`, err);
          } else {
            console.error(`[Backend] Failed to load adata.obsm[${columnName}][${index}]`, err);
          }
        },
      });
  }

  /**
   * Retrieve the computed CSS variable value for --ftable-min-width
   */
  getComputedMinWidth(): string {
    const computedStyle = getComputedStyle(this.el.nativeElement);
    const minWidth = computedStyle.getPropertyValue('--ftable-min-width').trim();
    return minWidth || '0';
  }


  // --- ChromVAR sum state ---
  private chromvarBaseMotif: string | null = null;
  private sumMotifs = new Set<string>();

  private getMotifId(row: any): string {
    return String(row?.motif_id ?? row?.['motif_id']);
  }

  isSumSelected(row: any): boolean {
    const id = this.getMotifId(row);
    return this.sumMotifs.has(id);
  }

  toggleSumMotif(row: any, checked: boolean): void {
    const id = this.getMotifId(row);
    if (!id) return;

    if (checked) this.sumMotifs.add(id);
    else this.sumMotifs.delete(id);

    // If user already clicked "Show" for a base motif, update the displayed score immediately
    this.updateChromvarCombinedIfReady();
  }

  async onShowAction(action: string, row: any): Promise<void> {
    const geneName = String(row.index);

    // DGEA-only UI action:
    // do NOT call backend here, let the parent convert it to gene_expression
    if (action === 'show_on_plot') {
      if (this.isCompare) {
        this.geneSelectedCompare.emit({ gene: geneName, action });
      } else {
        this.geneSelected.emit({ gene: geneName, action });
      }
      return;
    }

    if (action === 'chromvar_spot_scores') {
      this.chromvarBaseMotif = this.getMotifId(row);
      this.updateChromvarCombinedIfReady(true);

      if (this.isCompare) {
        this.geneSelectedCompare.emit({ gene: geneName, action });
      } else {
        this.geneSelected.emit({ gene: geneName, action });
      }
      return;
    }

    if (this.isCompare) {
      this.geneSelectedCompare.emit({ gene: geneName, action });
    } else {
      this.geneSelected.emit({ gene: geneName, action });
    }

    this.fetchAndUpdate(action, geneName);
  }

  private updateChromvarCombinedIfReady(force = false): void {
    if (!this.chromvarBaseMotif) {
      // Nothing to show yet until user clicks "Show"
      return;
    }

    // combined list = base motif + all checked motifs
    const motifs = new Set<string>([this.chromvarBaseMotif, ...this.sumMotifs]);

    // If you want "Show" to always hit backend even with no boxes ticked, keep force=true
    if (!force && motifs.size === 0) return;

    const index = Array.from(motifs).join(','); // backend will split and sum
    void this.fetchAndUpdate('chromvar_spot_scores', index);
  }

  get selectedSumCount(): number {
    return this.sumMotifs.size;
  }

  showSumChromvar(): void {
    const motifs = Array.from(this.sumMotifs);
    if (!motifs.length) return;

    //  sum-only (ignores base motif)
    void this.fetchAndUpdate('chromvar_spot_scores', motifs.join(','));

  }


}
