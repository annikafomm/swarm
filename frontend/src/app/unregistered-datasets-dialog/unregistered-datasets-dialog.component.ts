import { Component, OnInit, OnDestroy } from '@angular/core';
import { CommonModule } from '@angular/common';
import { MatDialogModule, MatDialogRef } from '@angular/material/dialog';
import { MatButtonModule } from '@angular/material/button';
import { MatIconModule } from '@angular/material/icon';
import { MatProgressSpinnerModule } from '@angular/material/progress-spinner';
import { MatTableModule } from '@angular/material/table';
import { MatCheckboxModule } from '@angular/material/checkbox';
import { MatCardModule } from '@angular/material/card';
import { MatTooltipModule } from '@angular/material/tooltip';
import { Subject } from 'rxjs';
import { takeUntil } from 'rxjs/operators';
import { DatasetService } from '../datasets.service';

interface UnregisteredDataset {
  dataset_id: string;
  alias: string;
  dataset_type: string;
  created_at: string;
  config_path: string;
  adata_path: string;
  adata_exists: boolean;
  status: 'ready' | 'missing_files';
  missing_files?: string | null;
  is_complete?: boolean;
}

@Component({
  selector: 'app-unregistered-datasets-dialog',
  standalone: true,
  imports: [
    CommonModule,
    MatDialogModule,
    MatButtonModule,
    MatIconModule,
    MatProgressSpinnerModule,
    MatTableModule,
    MatCheckboxModule,
    MatCardModule,
    MatTooltipModule,
  ],
  templateUrl: './unregistered-datasets-dialog.component.html',
  styleUrl: './unregistered-datasets-dialog.component.scss',
})
export class UnregisteredDatasetsDialogComponent implements OnInit, OnDestroy {
  datasets: UnregisteredDataset[] = [];
  isLoading = true;
  selectedDatasets = new Set<string>();
  actionInProgress = new Set<string>();
  private initialDatasetCount = 0; // Track initial count
  private destroy$ = new Subject<void>();

  displayedColumns: string[] = ['checkbox', 'alias', 'dataset_type', 'created_at', 'status', 'actions'];

  constructor(
    private datasetService: DatasetService,
    public dialogRef: MatDialogRef<UnregisteredDatasetsDialogComponent>
  ) { }

  ngOnInit(): void {
    this.loadUnregisteredDatasets();
  }

  ngOnDestroy(): void {
    this.destroy$.next();
    this.destroy$.complete();
  }

  private loadUnregisteredDatasets(): void {
    this.isLoading = true;
    this.datasetService.loadUnregisteredDatasets()
      .pipe(takeUntil(this.destroy$))
      .subscribe({
        next: (response) => {
          this.datasets = response.datasets || [];
          this.initialDatasetCount = this.datasets.length;
          this.isLoading = false;
        },
        error: (error) => {
          console.error('Failed to load unregistered datasets:', error);
          this.isLoading = false;
          this.datasets = [];
        }
      });
  }

  toggleDatasetSelection(datasetId: string): void {
    const dataset = this.datasets.find(d => d.dataset_id === datasetId);
    if (dataset && !dataset.is_complete) {
      return; // Don't allow selection of incomplete datasets
    }
    if (this.selectedDatasets.has(datasetId)) {
      this.selectedDatasets.delete(datasetId);
    } else {
      this.selectedDatasets.add(datasetId);
    }
  }

  isDatasetSelected(datasetId: string): boolean {
    return this.selectedDatasets.has(datasetId);
  }

  toggleAllSelection(): void {
    // Only allow selecting complete datasets
    const selectableDatasets = this.datasets.filter(d => d.is_complete);
    if (this.selectedDatasets.size === selectableDatasets.length) {
      this.selectedDatasets.clear();
    } else {
      selectableDatasets.forEach(d => this.selectedDatasets.add(d.dataset_id));
    }
  }

  isAllSelected(): boolean {
    const selectableDatasets = this.datasets.filter(d => d.is_complete);
    return selectableDatasets.length > 0 && this.selectedDatasets.size === selectableDatasets.length;
  }

  registerDataset(dataset: UnregisteredDataset): void {
    if (!dataset.is_complete) {
      alert(`Cannot register: Dataset is incomplete.\n\nMissing files:\n${dataset.missing_files || 'Unknown'}\n\nPlease ensure all required output files are present.`);
      return;
    }

    this.actionInProgress.add(dataset.dataset_id);
    this.datasetService.registerDataset(dataset.dataset_id)
      .pipe(takeUntil(this.destroy$))
      .subscribe({
        next: (response) => {
          console.log('✓ Dataset registered:', response);
          // Remove from list
          this.datasets = this.datasets.filter(d => d.dataset_id !== dataset.dataset_id);
          this.selectedDatasets.delete(dataset.dataset_id);
          this.actionInProgress.delete(dataset.dataset_id);
        },
        error: (error) => {
          console.error('Failed to register dataset:', error);
          alert(`Error: ${error.error?.detail || 'Failed to register dataset'}`);
          this.actionInProgress.delete(dataset.dataset_id);
        }
      });
  }

  deleteDataset(dataset: UnregisteredDataset): void {
    if (!confirm(`Delete dataset "${dataset.alias}" and its files?`)) {
      return;
    }

    this.actionInProgress.add(dataset.dataset_id);
    this.datasetService.deleteUnregisteredDataset(dataset.dataset_id)
      .pipe(takeUntil(this.destroy$))
      .subscribe({
        next: (response) => {
          console.log('✓ Dataset deleted:', response);
          // Remove from list
          this.datasets = this.datasets.filter(d => d.dataset_id !== dataset.dataset_id);
          this.selectedDatasets.delete(dataset.dataset_id);
          this.actionInProgress.delete(dataset.dataset_id);
        },
        error: (error) => {
          console.error('Failed to delete dataset:', error);
          alert(`Error: ${error.error?.detail || 'Failed to delete dataset'}`);
          this.actionInProgress.delete(dataset.dataset_id);
        }
      });
  }

  registerSelectedDatasets(): void {
    if (this.selectedDatasets.size === 0) {
      alert('Please select at least one dataset to register');
      return;
    }

    const selectedArray = Array.from(this.selectedDatasets);
    let completedCount = 0;

    selectedArray.forEach(datasetId => {
      const dataset = this.datasets.find(d => d.dataset_id === datasetId);
      if (dataset) {
        this.registerDataset(dataset);
        // Count completion
        this.actionInProgress.add(datasetId);
      }

      // Check if all completed
      setTimeout(() => {
        completedCount++;
        if (completedCount === selectedArray.length) {
          console.log('All registrations completed');
        }
      }, 100);
    });
  }

  deleteSelectedDatasets(): void {
    if (this.selectedDatasets.size === 0) {
      alert('Please select at least one dataset to delete');
      return;
    }

    if (!confirm(`Delete ${this.selectedDatasets.size} dataset(s) and their files?`)) {
      return;
    }

    const selectedArray = Array.from(this.selectedDatasets);
    selectedArray.forEach(datasetId => {
      const dataset = this.datasets.find(d => d.dataset_id === datasetId);
      if (dataset) {
        this.deleteDataset(dataset);
      }
    });
  }

  close(): void {
    const datasetsChanged = this.datasets.length < this.initialDatasetCount;
    this.dialogRef.close({ datasetsChanged });
  }

  getStatusBadgeClass(status: string): string {
    return status === 'ready' ? 'status-badge status-ready' : 'status-badge status-missing';
  }

  getStatusTooltip(dataset: UnregisteredDataset): string {
    if (dataset.is_complete) {
      return 'All required files present. Ready to register.';
    }
    return `Incomplete: Missing ${dataset.missing_files || 'required files'}`;
  }

  isDatasetIncomplete(dataset: UnregisteredDataset): boolean {
    return !dataset.is_complete;
  }
}
