import { Component } from '@angular/core';
import { FormBuilder, FormGroup, ReactiveFormsModule, FormControl } from '@angular/forms';
import { CommonModule } from '@angular/common';
import { Router } from '@angular/router'

@Component({
  selector: 'app-form-page',
  standalone: true,
  imports: [ReactiveFormsModule, CommonModule],
  templateUrl: './form-page.component.html',
  styleUrls: ['./form-page.component.scss'],
})
export class FormPageComponent {
  form: FormGroup;
  uploadedFile: File | null = null;
  singleCellFile: File | null = null;
  uploadedFile2: File | null = null;
  spongeNetworkFile: File | null = null;
  genieFile: File | null = null;

  availableScores = ['liana+', 'spongeeffects', 'squidpy', 'viper', 'aucell'];

  constructor(private fb: FormBuilder, private router: Router) {
    this.form = this.fb.group({
      dataset: ['Visium'],
      method: ['Genie3'],
      normalization: [false],
      filteringSpatial: [false],
      filteringSingleCell: [false],
      tangram: [false],
      scores: this.fb.group({
        'liana+': [false],
        'spongeeffects': [false],
        'squidpy': [false],
        'viper': [false],
        'aucell': [false],
      })
    });
  }

  // Upload-Funktionen
  onFileSelected(event: Event, fileType: 'spatial' | 'singleCell' | 'file2' | 'sponge' | 'genie') {
    const input = event.target as HTMLInputElement;
    if (!input.files || input.files.length === 0) return;

    const file = input.files[0];

    switch (fileType) {
      case 'spatial':
        if (!file.name.endsWith('.h5ad')) return;
        this.uploadedFile = file;
        break;
      case 'singleCell':
        if (!file.name.endsWith('.h5ad')) return;
        this.singleCellFile = file;
        break;
      case 'file2':
        this.uploadedFile2 = file;
        break;
      case 'sponge':
        this.spongeNetworkFile = file;
        break;
      case 'genie':
        this.genieFile = file;
        break;
    }
  }

  // Check, ob eine Score-Abhängigkeit vorliegt
  requiresSponge(): boolean {
    const scores = this.form.get('scores')?.value;
    const spongeNeeded = scores['spongeeffects'] || scores['aucell'];
    // Wenn file2 (oben) schon hochgeladen wurde → Sponge nicht mehr erforderlich
    return spongeNeeded && !this.uploadedFile2;
  }

  requiresGenie(): boolean {
    const scores = this.form.get('scores')?.value;
    const genieNeeded = scores['viper'] || scores['aucell'];
    // Wenn file2 (oben) schon hochgeladen wurde → Genie nicht mehr erforderlich
    return genieNeeded && !this.uploadedFile2;
  }


  onSubmit() {
    const formData = new FormData();

    formData.append('dataset', this.form.value.dataset);
    formData.append('method', this.form.value.method);
    formData.append('normalization', this.form.value.normalization.toString());
    formData.append('filteringSpatial', this.form.value.filteringSpatial.toString());
    formData.append('filteringSingleCell', this.form.value.filteringSingleCell.toString());
    formData.append('tangram', this.form.value.tangram.toString());

    // Scores als JSON
    formData.append('scores', JSON.stringify(this.form.value.scores));

    if (this.uploadedFile) formData.append('spatialFile', this.uploadedFile);
    if (this.singleCellFile) formData.append('singleCellFile', this.singleCellFile);
    if (this.uploadedFile2) formData.append('file2', this.uploadedFile2);
    if (this.spongeNetworkFile) formData.append('spongeNetwork', this.spongeNetworkFile);
    if (this.genieFile) formData.append('genieFile', this.genieFile);

    console.log('Form Data:', this.form.value);
    this.router.navigate(['/plot'])
  }
}
