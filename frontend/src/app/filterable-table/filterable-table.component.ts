import { Component, OnInit, Input, Output, EventEmitter, OnChanges, SimpleChanges } from '@angular/core';
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
  @Output() featuresUpdated = new EventEmitter<void>();

  constructor(
    private http: HttpClient,
    private sessionService: SessionService,
  ) { }

  columns: string[] = [];
  rows: any[] = [];
  filters: { [col: string]: string } = {};
  sortColumn: string | null = null;
  sortAsc: boolean = true;

  // Pagination
  pageSize = 50; // number of rows per page
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
    if (Array.isArray(this.data)) {
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

  async fetchAndUpdate(columnName: string, index: string) {
    let isGeneExpression = columnName === 'gene_expression';
    let request = isGeneExpression
      ? `${this.sessionService.apiUrl}/X/${index}`
      : `${this.sessionService.apiUrl}/obsm/${columnName}/${index}`;
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
          } else {
            console.log(`[Backend] Loaded adata.obsm[${columnName}][${index}]`);
          }

          this.featuresUpdated.emit();
        },
        error: (err) => {
          if (isGeneExpression) {
            console.error(`[Backend] Failed to load adata[:, ${index}].X`, err);
          } else {
            console.error(
              `[Backend] Failed to load adata.obsm[${columnName}][${index}]`,
              err,
            );
          }
        },
      });
  }
}
