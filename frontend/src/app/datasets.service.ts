import { Injectable, Injector } from '@angular/core';
import { BehaviorSubject, Observable } from 'rxjs';
import { HttpClient } from '@angular/common/http';
import { SessionService } from './session.service';

export interface Dataset {
  id: string;
  alias: string;
  adata_path: string;
  tangram_adata_path?: string;
  genie_network_path?: string;
  sponge_network_path?: string;
  geojson_path?: string;
  type: 'builtin' | 'uploaded';
  created_at?: string;
}

@Injectable({
  providedIn: 'root'
})
export class DatasetService {
  private availableDatasets = new BehaviorSubject<Dataset[]>([]);
  public availableDatasets$ = this.availableDatasets.asObservable();

  // Main dataset for visualization
  private selectedDataset = new BehaviorSubject<Dataset | null>(null);
  public selectedDataset$ = this.selectedDataset.asObservable();

  // Comparison dataset (for side-by-side analysis)
  private selectedDatasetCompare = new BehaviorSubject<Dataset | null>(null);
  public selectedDatasetCompare$ = this.selectedDatasetCompare.asObservable();

  constructor(private http: HttpClient, private injector: Injector) {
    // Datasets will be loaded after session is initialized
  }

  loadAvailableDatasets(): void {
    this.http.get<any>('/api/datasets', { withCredentials: true }).subscribe({
      next: (response) => {
        const datasets = this.flattenDatasets(response);
        this.availableDatasets.next(datasets);
      },
      error: (err) => console.error('Failed to load datasets', err)
    });
  }

  private flattenDatasets(registry: any): Dataset[] {
    const all: Dataset[] = [];

    // Add builtin datasets
    if (registry.builtin) {
      Object.entries(registry.builtin).forEach(([id, data]: [string, any]) => {
        all.push({ id, ...data, type: 'builtin' });
      });
    }

    // Add uploaded datasets
    if (registry.uploaded) {
      Object.entries(registry.uploaded).forEach(([id, data]: [string, any]) => {
        all.push({ id, ...data, type: 'uploaded' });
      });
    }

    return all;
  }

  // Main dataset selection
  selectDataset(dataset: Dataset | null): void {
    this.selectedDataset.next(dataset);
    if (dataset) {
      this.loadNetworksForDataset(dataset);
    }
  }

  getSelectedDataset(): Dataset | null {
    return this.selectedDataset.value;
  }

  // Comparison dataset selection
  selectDatasetCompare(dataset: Dataset | null): void {
    this.selectedDatasetCompare.next(dataset);
    if (dataset) {
      this.loadNetworksForDataset(dataset);
    }
  }

  getSelectedDatasetCompare(): Dataset | null {
    return this.selectedDatasetCompare.value;
  }

  private loadNetworksForDataset(dataset: Dataset): void {
    // Lazy-inject SessionService to avoid circular dependency
    const sessionService = this.injector.get(SessionService);

    if (dataset.genie_network_path) {
      this.http.post(
        `${sessionService.apiUrl}/read_network_genie`,
        { path: dataset.genie_network_path },
        { withCredentials: true }
      ).subscribe({
        next: () => console.log('✓ Loaded Genie3 network'),
        error: (err) => console.error('✗ Failed to load Genie3 network:', err)
      });
    }

    if (dataset.sponge_network_path) {
      this.http.post(
        `${sessionService.apiUrl}/read_network_sponge`,
        { path: dataset.sponge_network_path },
        { withCredentials: true }
      ).subscribe({
        next: () => console.log('✓ Loaded SPONGE network'),
        error: (err) => console.error('✗ Failed to load SPONGE network:', err)
      });
    }
  }
}
