// src/app/form-page/form-page.component.ts

// import for an Angular component
import { Component } from '@angular/core';
// Imports the CommonModule (provides core Angular directives like *ngIf, *ngFor).
import { CommonModule } from '@angular/common';
/* Imports classes for Reactive Forms:
   FormBuilder: helper to create form groups.
   FormGroup: represents a collection of form controls.
   ReactiveFormsModule: enables Reactive Forms in Angular. */
import { FormBuilder, FormGroup, ReactiveFormsModule, Validators } from '@angular/forms';
import { Router } from '@angular/router';
// Imports Angular's HTTP client for API requests
import { HttpClient, HttpEventType, HttpClientModule } from '@angular/common/http';
// Imports environment variables, e.g., the backend API URL.
import { environment } from '../../environments/environment';

/* Defines the Angular component:
   selector: Tag name in HTML (<app-form-page>).
   standalone: true — the component is standalone and does not need to be declared in a module.
   imports: modules this component depends on.
   templateUrl: path to the HTML template. */
@Component({
  selector: 'app-form-page',
  standalone: true,
  imports: [CommonModule, ReactiveFormsModule, HttpClientModule],
  templateUrl: './form-page.component.html',
  styleUrls: ['./form-page.component.scss'],
})

// Begins the component class definition
export class FormPageComponent {
  // Form object
  form: FormGroup;

  // Files
  spatialFile!: File;        // the spatial file is required
  singleCellFile?: File;
  precomputedFile?: File;
  spongeNetworkFile?: File;
  genieFile?: File;

  // This is relevant for the HTML view — which scores should be displayed/available
  availableScores = ['LIANA+', 'SPONGeffects', 'squidpy', 'VIPER', 'AUCell'];

  uploading = false;                  // Upload is running; set to true as soon as submit starts
  uploadProgress: number | null = null; // used for a potential upload progress bar
  errorMsg: string | null = null;       // used for error messages

  constructor(
    private fb: FormBuilder,   // to easily build reactive forms
    private router: Router,    // important to navigate to another page later
    private http: HttpClient   // to send API requests to the backend
  ) {
    // Build the reactive form
    this.form = this.fb.group({
      dataset: ['Visium', Validators.required],           // dataset name or identifier
      method: ['Genie3', Validators.required],            // selected method (e.g., spatial)
      tangram: [false],                             // whether Tangram alignment is activated

      normalization: [false],
      filteringSpatial: [false],
      filteringSingleCell: [false],

      // group of score checkboxes; the key is the score name, the value is a boolean
      scores: this.fb.group({
        'LIANA+': [false],
        'SPONGeffects': [false],
        'squidpy': [false],
        'VIPER': [false],
        'AUCell': [false],
      }),
    });
  }

  // Convenience getter for the score checkbox group value
  private get scoreValues(): Record<string, boolean> {
    return (this.form.get('scores')?.value as Record<string, boolean>) || {};
  }

  // Returns the names of all selected scores
  private get selectedScores(): string[] {
    const s = this.scoreValues;
    return Object.keys(s).filter((k) => !!s[k]);
  }

  // File selection handlers --------------------------------------------------

  // Spatial dataset file (required)
  onSpatialSelected(event: Event): void {
    const file = (event.target as HTMLInputElement)?.files?.[0] || null;
    if (file) this.spatialFile = file;
  }

  // Single-cell reference file (required if Tangram is enabled)
  onSingleCellSelected(event: Event): void {
    const file = (event.target as HTMLInputElement)?.files?.[0] || null;
    this.singleCellFile = file ?? undefined;
  }

  // Precomputed file (if provided, network and genie files may not be required)
  onPrecomputedSelected(event: Event): void {
    const file = (event.target as HTMLInputElement)?.files?.[0] || null;
    this.precomputedFile = file ?? undefined;
  }

  // SPONGE network file — required if SPONGeffects or AUCell are selected (and no precomputed file is provided)
  onSpongeNetworkSelected(event: Event): void {
    const file = (event.target as HTMLInputElement)?.files?.[0] || null;
    this.spongeNetworkFile = file ?? undefined;
  }

  // GENIE/VIPER network file — required if VIPER or AUCell are selected (and no precomputed file is provided)
  onGenieSelected(event: Event): void {
    const file = (event.target as HTMLInputElement)?.files?.[0] || null;
    this.genieFile = file ?? undefined;
  }

  onFileSelected(event: Event, type: string): void {
    const input = event.target as HTMLInputElement;
    if (input?.files && input.files.length > 0) {
      const file = input.files[0];

      switch (type) {
        case 'spatial':
          this.spatialFile = file;
          break;
        case 'singleCell':
          this.singleCellFile = file;
          break;
        case 'precomputed':
          this.precomputedFile = file;
          break;
        case 'sponge':
          this.spongeNetworkFile = file;
          break;
        case 'genie':
          this.genieFile = file;
          break;
      }
    }
  }


  // Toggle a score programmatically (e.g., from the template)
  changeScore(scoreName: string, checked: boolean): void {
    const ctrl = (this.form.get('scores') as FormGroup)?.get(scoreName);
    if (ctrl) ctrl.setValue(checked);
  }

  // Whether Tangram is required
  requiresTangram(): boolean {
    return !!this.form.value.tangram && !this.singleCellFile; // if tangram is enabled, a single-cell file must be present
  }

  // Whether a SPONGE network file is required
  requiresSponge(): boolean {
    const s = this.scoreValues; // gets the current value of the scores checkbox group
    const spongeNeeded = s['SPONGeffects'] || s['AUCell']; // true if one of these two is selected
    return spongeNeeded && !this.precomputedFile; // but only if no precomputed file was uploaded
  }

  // Same logic for the GENIE/VIPER file
  requiresGenie(): boolean {
    const s = this.scoreValues;
    const needsGenie = s['VIPER'] || s['AUCell'] || s['SPONGeffects'];
    return needsGenie && !this.precomputedFile;
  }


  // Final client-side validation before enabling submit
  canSubmit(): boolean {
    if (!this.spatialFile) return false; // is there a spatial dataset?
    if (this.form.value.tangram && !this.singleCellFile) return false; // if tangram is enabled, a single-cell file must be uploaded
    if (this.requiresSponge() && !this.spongeNetworkFile) return false; // same as for Tangram
    if (this.requiresGenie() && !this.genieFile) return false;          // same logic
    return this.form.valid;
  }

  shouldOfferGenieUpload(): boolean {
    const s = this.scoreValues;
    return !!s['LIANA+'];
  }

  onSubmit(): void {
    // 1) Final client-side validation
    if (!this.canSubmit()) {
      // prevent incomplete data from being sent forward
      return this.fail('Please complete all required fields/files before submitting.');
    }

    // 2) Build FormData only after validation
    // here the multipart/form-data payload is created
    const fd = new FormData();
    fd.append('dataset', this.form.value.dataset);
    fd.append('method', this.form.value.method);
    fd.append('tangram', String(!!this.form.value.tangram));

    fd.append('normalization', String(!!this.form.value.normalization));
    fd.append('filteringSpatial', String(!!this.form.value.filteringSpatial));
    fd.append('filteringSingleCell', String(!!this.form.value.filteringSingleCell));

    // Append scores as JSON array of selected names
    fd.append('scores', JSON.stringify(this.selectedScores));

    // Required file
    fd.append('spatialFile', this.spatialFile);

    // Optional files depending on toggles/selection
    if (this.singleCellFile) fd.append('singleCell', this.singleCellFile);
    if (this.precomputedFile) fd.append('precomputed', this.precomputedFile);
    if (this.spongeNetworkFile) fd.append('spongeNetwork', this.spongeNetworkFile);
    


    if (this.genieFile) {
      fd.append('genieNetwork', this.genieFile);
    }

    const s = this.scoreValues;
  if (s['LIANA+']) {
    // Wenn keine Datei hochgeladen wurde, soll Backend Default nehmen
    fd.append('useDefaultLiana', this.genieFile ? 'false' : 'true');
  }



    // 3) Send the request
    this.uploading = true;
    this.uploadProgress = 0;

    this.http
      .post(`${environment.apiBaseUrl}/api/upload`, fd, {
        reportProgress: true, // so we can show upload progress later
        observe: 'events',
      })
      .subscribe({
        next: (event) => {
          // Update the progress bar (if any)
          if (event.type === HttpEventType.UploadProgress) {
            const total = (event.total ?? 0) || 0;
            // Guard against division by zero
            this.uploadProgress = total ? Math.round((100 * (event.loaded ?? 0)) / total) : null;
          }

          // When the response arrives
          if (event.type === HttpEventType.Response) {
            this.uploading = false;
            this.uploadProgress = null;

            // Navigate to the results page with the returned jobId
            // (expects the backend to respond with { jobId: string })
            const jobId = (event.body as any)?.jobId;
            if (!jobId) {
              return this.fail('Upload succeeded, but no jobId was returned by the server.');
            }
            this.router.navigate(['/result', jobId]);
          }
        },
        error: (err) => {
          this.uploading = false;
          this.uploadProgress = null;
          // Übergib das ganze err.error – normalizeError macht den Rest
          this.fail(err?.error ?? err ?? 'Upload failed.');
        },
      });
  }

  private normalizeError(e: any): string {
    if (!e) return 'Unknown error';
    if (typeof e === 'string') return e;
    if (Array.isArray(e)) {
      // Mappe Arrays (z. B. [{msg:'a'},{msg:'b'}]) auf eine lesbare Liste
      return e.map(item => this.normalizeError(item)).join('\n');
    }
    if (typeof e.message === 'string') return e.message;
    if (typeof e.detail === 'string') return e.detail;
    try { return JSON.stringify(e); } catch { return String(e); }
  }


  // For showing error messages with auto-clear
  private fail(msg: any): void {
    this.errorMsg = this.normalizeError(msg);
    setTimeout(() => (this.errorMsg = null), 10000);
  }

}
