import { Component, OnInit, Input, Output, EventEmitter } from '@angular/core';
import { FormsModule } from '@angular/forms';
import { CommonModule } from '@angular/common';
import { SessionService } from '../session.service';
import { HttpClient } from '@angular/common/http';

@Component({
  selector: 'app-table',
  templateUrl: './filterable-table.component.html',
  styleUrls: ['./filterable-table.component.scss'],
  standalone: true,
  imports: [FormsModule, CommonModule],
})
export class FilterableTableComponent implements OnInit {
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
  ) {}

  columns: string[] = [];
  rows: any[] = [];
  filters: { [col: string]: string } = {};
  sortColumn: string | null = null;
  sortAsc: boolean = true;

  ngOnInit() {
    this.prepareTable();
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

  displayNumeric(value: any): string {
    if (typeof value === 'number') {
      return value.toFixed(3); // display 3 decimal places, could change
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
    this.sessionService
      .callWithSession(() =>
        this.http.get(
          `${this.sessionService.apiUrl}/obsm/${columnName}/${index}`,
          { withCredentials: true },
        ),
      )
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
          console.log(`[Backend] Loaded adata.obsm["${columnName}][${index}]`);
          this.featuresUpdated.emit();
        },
        error: (err) =>
          console.error(
            `[Backend] Failed to load adata.obsm["${columnName}][${index}]`,
            err,
          ),
      });
  }
}
