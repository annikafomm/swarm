import { Component, OnDestroy } from '@angular/core';
import { Location } from '@angular/common';
import { CommonModule } from '@angular/common';
import { Subject } from 'rxjs';
import { takeUntil } from 'rxjs/operators';
import { DatasetService } from '../datasets.service';
import { Dataset } from '../datasets.service';

@Component({
  selector: 'app-info',
  standalone: true,
  imports: [CommonModule],
  templateUrl: './info.component.html',
  styleUrl: './info.component.scss'
})
export class InfoComponent implements OnDestroy {
  selectedDataset: Dataset | null = null;
  private destroy$ = new Subject<void>();

  constructor(
    private location: Location,
    private datasetService: DatasetService
  ) {
    // Subscribe to dataset selection changes
    this.datasetService.selectedDataset$
      .pipe(takeUntil(this.destroy$))
      .subscribe(dataset => {
        this.selectedDataset = dataset;
      });
  }

  goBack() {
    this.location.back();
  }

  /**
   * Check if any multiome features are available in the selected dataset
   */
  hasMultiomeFeatures(): boolean {
    if (!this.selectedDataset) return false;
    const ds = this.selectedDataset;
    return !!(
      ds.use_chromvar ||
      ds.use_differential_motif_activity ||
      ds.use_footprinting ||
      ds.use_moranI ||
      ds.use_gearyC
    );
  }

  /**
   * Get data file paths as an array for iteration in template
   */
  getMultiomeFiles(): Array<{ name: string; path: string; label: string }> {
    if (!this.selectedDataset) return [];
    const ds = this.selectedDataset;
    const files: Array<{ name: string; path: string | undefined; label: string }> = [
      { name: 'ST Scores', path: ds.adata_st_scores_path, label: 'adata_st_scores' },
      { name: 'TG Scores (Tangram)', path: ds.adata_tg_scores_path, label: 'adata_tg_scores' },
      { name: 'Tangram Mapping', path: ds.adata_map_path, label: 'adata_map' },
      { name: 'Tangram Mapping (CSV)', path: ds.adata_map_X_csv_path, label: 'adata_map.X' },
      { name: 'Spot Metadata (CSV)', path: ds.adata_map_var_csv_path, label: 'adata_map.var' },
      { name: 'Pipeline Log', path: ds.calc_scores_log_path, label: 'calc_scores' },
      { name: 'Seurat Multiome Object', path: ds.global_motif_analysis_path, label: 'global_motif_analysis' },
      { name: 'Motif-to-TF Mapping', path: ds.motif_to_tf_csv_path, label: 'motif_to_tf' },
      { name: 'ChromVAR Scores', path: ds.chromvar_scores_csv_path, label: 'chromvar_scores' },
    ];
    return files.filter((f): f is { name: string; path: string; label: string } => !!f.path);
  }

  /**
   * Get downloadable CSVs from differential motif activity analysis
   */
  getDifferentialMotifActivityFiles(): Array<{ comparison: string; path: string }> {
    if (!this.selectedDataset?.diff_motif_activity_csv_paths) return [];
    return Object.entries(this.selectedDataset.diff_motif_activity_csv_paths).map(
      ([comparison, path]) => ({ comparison, path: path as string })
    );
  }

  /**
   * Get downloadable footprint PDFs
   */
  getFootprintPDFs(): Array<{ motifId: string; path: string }> {
    if (!this.selectedDataset?.footprint_pdf_paths) return [];
    return Object.entries(this.selectedDataset.footprint_pdf_paths).map(
      ([motifId, path]) => ({ motifId, path: path as string })
    );
  }

  /**
   * Generate download link for a file
   */
  getDownloadLink(filePath: string): string {
    // Encode the path for URL parameter
    return `/api/download/${encodeURIComponent(filePath)}`;
  }

  /**
   * Format creation date to readable string
   */
  formatDate(dateStr?: string): string {
    if (!dateStr) return 'Unknown';
    try {
      return new Date(dateStr).toLocaleDateString();
    } catch {
      return dateStr;
    }
  }

  ngOnDestroy() {
    this.destroy$.next();
    this.destroy$.complete();
  }

}
