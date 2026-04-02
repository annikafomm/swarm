import { Component, OnInit, Input, Output, EventEmitter, OnChanges, SimpleChanges, ElementRef } from '@angular/core';
import { FormsModule } from '@angular/forms';
import { CommonModule } from '@angular/common';
import { SessionService } from '../session.service';
import { HttpClient } from '@angular/common/http';
import { TranslatePipe } from '../translate.pipe';

@Component({
  selector: 'app-table',
  templateUrl: './filterable-table.component.html',
  styleUrls: ['./filterable-table.component.scss'],
  standalone: true,
  imports: [FormsModule, CommonModule, TranslatePipe],
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
  @Output() featuresUpdated = new EventEmitter<void>();
  @Output() geneSelected = new EventEmitter<{ gene: string; action: string }>();
  @Output() geneSelectedCompare = new EventEmitter<{ gene: string; action: string }>();
  constructor(
    private http: HttpClient,
    private sessionService: SessionService,
    private el: ElementRef,
  ) { }

  columns: string[] = [];
  rows: any[] = [];
  filters: { [col: string]: string } = {};
  sortColumn: string | null = null;
  sortAsc: boolean = true;
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
      this.sortColumn = null;
      this.sortAsc = true;
      this.filters = {};
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
    if (this.sortColumn) {
      const sortKey = this.sortColumn as string;
      result.sort((a, b) => {
        const valA = a[sortKey];
        const valB = b[sortKey];

        if (typeof valA === 'number' && typeof valB === 'number') {
          return this.sortAsc ? valA - valB : valB - valA;
        }
        return this.sortAsc
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
    if (this.sortColumn === col) {
      this.sortAsc = !this.sortAsc;
    } else {
      this.sortColumn = col;
      this.sortAsc = true;
    }
  }

  setSort(col: string, ascending: boolean) {
    this.sortColumn = col;
    this.sortAsc = ascending;
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

          this.featuresUpdated.emit();
        },
        error: (err) => {
          if (this.features) {
            for (const feature of this.features) {
              feature.properties[this.updateColumn] = undefined;
            }
            this.featuresUpdated.emit();
          }

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
        this.geneSelectedCompare.emit({ gene: geneName, action: action });
      } else {
        this.geneSelected.emit({ gene: geneName, action: action });
      }
      return;
    }

    // Emit the gene selection event for all score column actions
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
