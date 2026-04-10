import { Injectable, Injector } from '@angular/core';
import { BehaviorSubject, Observable } from 'rxjs';
import { HttpClient } from '@angular/common/http';
import { SessionService } from './session.service';

export interface Dataset {
  // Basic fields
  id: string;
  alias: string;
  adata_path: string;
  type: 'builtin' | 'uploaded';
  user?: string;
  created_at?: string;
  dataset_type?: string; // e.g., "Visium", "Xenium", "Multiome"

  // Core spatial/network paths
  tangram_adata_path?: string;
  genie_network_path?: string;
  sponge_network_path?: string;
  geojson_path?: string;

  // Multiome: Always-present outputs
  adata_st_scores_path?: string;      // Spatial transcriptomics with scores (H5AD)
  adata_tg_scores_path?: string;      // Tangram-projected data with scores (H5AD)
  adata_map_path?: string;            // Tangram mapping matrix (cells × spots, H5AD)
  adata_map_X_csv_path?: string;      // Tangram mapping as CSV
  adata_map_var_csv_path?: string;    // Spot metadata from Tangram
  calc_scores_log_path?: string;      // Pipeline execution log
  global_motif_analysis_path?: string; // Complete Seurat multiome object (RDS)

  // Multiome: Conditional outputs
  motif_to_tf_csv_path?: string;      // Motif ID ↔ TF name mapping
  spot_obj_chromvar_path?: string;    // Spot-level chromVAR object (RDS)
  spot_obj_footprints_path?: string;  // Spot-level footprints (RDS)
  dissociated_obj_footprints_path?: string; // Dissociated cell footprints (RDS)
  chromvar_scores_csv_path?: string;  // chromVAR deviation scores (cells × motifs)
  diff_motif_activity_csv_paths?: { [comparison: string]: string }; // comparison → CSV path mapping
  footprint_pdf_paths?: { [motif_id: string]: string }; // motif_id → PDF path mapping

  // Multiome: Feature flags indicating which pipeline steps were executed
  use_chromvar?: boolean;
  use_differential_motif_activity?: boolean;
  use_footprinting?: boolean;
  use_moranI?: boolean;
  use_gearyC?: boolean;

  // Visium-specific flag
  use_tangram?: boolean;
  use_multiome?: boolean;

  // Deprecated: kept for backward compatibility
  footprint_list?: string[];
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

  // Load unregistered datasets from uploads folder
  loadUnregisteredDatasets(): Observable<any> {
    return this.http.get<any>('/api/unregistered_datasets', { withCredentials: true });
  }

  // Register an unregistered dataset
  registerDataset(datasetId: string): Observable<any> {
    const params = new URLSearchParams();
    params.set('dataset_id', datasetId);
    return this.http.post<any>(`/api/register_dataset?${params.toString()}`, {}, { withCredentials: true });
  }

  // Delete an unregistered dataset
  deleteUnregisteredDataset(datasetId: string): Observable<any> {
    const params = new URLSearchParams();
    params.set('dataset_id', datasetId);
    return this.http.post<any>(`/api/delete_unregistered_dataset?${params.toString()}`, {}, { withCredentials: true });
  }
}
