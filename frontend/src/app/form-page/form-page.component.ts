import { Component } from '@angular/core';
import { FormBuilder, FormGroup, Validators } from '@angular/forms';
import { CommonModule } from '@angular/common';
import { ReactiveFormsModule } from '@angular/forms';
import { HttpClientModule } from '@angular/common/http';
import { HttpClient, HttpEvent, HttpEventType } from '@angular/common/http';

import { PathsService } from '../paths.service';
import { DEFAULT_PATHS } from '../constants';

@Component({
  selector: 'app-form-page',
  standalone: true,
  imports: [CommonModule, ReactiveFormsModule, HttpClientModule],
  templateUrl: './form-page.component.html',
  styleUrls: ['./form-page.component.scss'],
})
export class FormPageComponent {

  form: FormGroup;

  // files
  spatialFile?: File;
  singleCellFile?: File;
  multiomeFile?: File;
  // network scores uploads (shared)
  genie3NetFile?: File;          // for VIPER and/or AUCell/GSVA/ssGSEA
  spongeNAFile?: File;           // SPONGE networkanalysis
  spongeNIFile?: File;           // SPONGE networkinteractions
  // LIANA+
  // LIANA+
  lianaGenie3File?: File;
  lianaPathwayFile?: File;


  // ui state
  uploading = false;
  uploadProgress = 0;
  errorMsg = '';

  resultJsonUrl: string | null = null;   // vom Server: Download-URL zur erzeugten JSON
  serverPayload: any = null;             // vom Server zurückgegebenes Payload-Objekt


  constructor(private fb: FormBuilder, private http: HttpClient, private pathsService: PathsService) {
    this.form = this.fb.group({

      email: [''],
      dataset: ['Visium', Validators.required],

      spatialOptions: this.fb.group({
        normalization: [false],
        filtering: [false],
      }),

      useTangram: [false],
      tangram: this.fb.group({
        filterSingleCell: [false],
        normalizeSingleCell: [false],
      }),

      useTangramMultiome: [false],
      // tangram: this.fb.group({
      //   filterSingleCell: [false],
      //   normalizeSingleCell: [false],
      // }),

      // score toggles
      scores: this.fb.group({
        networkScores: [false],
        squidpy: [false],
        lianaPlus: [false],
        chromVAR: [false],
        diffMotifActivity: [false],
        motifEnrichment: [false],
        FootprintingBias: [false]
      }),
      genome: ['hg38'],

      // network scores options
      network: this.fb.group({
        algorithms: this.fb.group({
          viper: [false],
          aucell: [false],
          gsva: [false],
          ssgsea: [false],
        }),
      params: this.fb.group({
        mScoreThreshold: ['0.1'],
        pAdjust: ['0.05'],
        ensemblId: ['ensembl_id'],
        featureCol: ['feature_type'],
        rnaTypes: ['lncRNA,protein_coding'],
        maxModules: ['20'],
        }),
      genie3Params: this.fb.group({
        topNWeights: [100000],
        nRegulatoryGenes: [20],
        nRegulons: [20],
        }),
      }),

      // squidpy options (fields are reused depending on method)
      squidpy: this.fb.group({
        methods: this.fb.group({
          moranI: [false],
          gearyC: [false],
          centrality_score: [false],
          co_occurrence: [false],
          neighborhood_enrichment: [false],
        }),
        moranI: this.fb.group({
          nPerms: [''],
          twoTailed: [false],
          corrMethod: ['fdr_bh'],
        }),
        gearyC: this.fb.group({
          nPerms: [''],
          twoTailed: [false],
          corrMethod: ['fdr_bh'],
        }),
        centrality_score: this.fb.group({
          clusterKey: ['leiden'],
        }),
        co_occurrence: this.fb.group({
          clusterKey: ['leiden'],
          interval: ['50'],
          nSplits: [''],
        }),
        neighborhood_enrichment: this.fb.group({
          clusterKey: ['leiden'],
          libraryKey: [''],
          nPerms: ['1000'],
        }),
      }),

      // LIANA+
      liana: this.fb.group({
        compositionColumn: ['tangram_ct_pred'],
      }),
    });

    // clear single-cell if Tangram toggled off
    this.form.get('useTangram')!.valueChanges.subscribe(on => {
      if (!on) {
        this.singleCellFile = undefined;
        this.form.get('tangram')!.reset({
          filterSingleCell: false,
          normalizeSingleCell: false,
        });
      }
    });

    // if network scores toggled off → reset inner state
    this.form.get('scores.networkScores')!.valueChanges.subscribe(on => {
      if (!on) {
        this.form.get('network.algorithms')!.reset({
          viper: false, aucell: false, gsva: false, ssgsea: false
        });
        this.genie3NetFile = undefined;
        this.spongeNAFile = undefined;
        this.spongeNIFile = undefined;
      }
    });
  }

  // ---------- getters for template type-safety ----------
  get tangramGroup() {
    return this.form.get('tangram') as FormGroup;
  }
  get squidpyGroup() {
    return this.form.get('squidpy') as FormGroup;
  }
  get usingSponge(): boolean {
    return !!(this.spongeNAFile && this.spongeNIFile);
  }
  get genie3Selected(): boolean {
    return !!this.genie3NetFile;
  }
  useTangramChecked(): boolean {
    return !!this.form.get('useTangram')?.value;
  }
  useTangramMultiomeChecked(): boolean {
    return !!this.form.get('useTangramMultiome')?.value;
  }
  squidpyMethodIs(m: string): boolean {
    return this.form.get('squidpy.method')?.value === m;
  }
  hasTriadSelected(): boolean {
    const a = this.form.get('network.algorithms.aucell')?.value;
    const g = this.form.get('network.algorithms.gsva')?.value;
    const s = this.form.get('network.algorithms.ssgsea')?.value;
    return !!(a || g || s);
  }

  // ---------- file handling ----------
  onFileSelected(evt: Event, type:
      'spatial' | 'singleCell' | 'multiome' | 'genie3Net' | 'spongeNA' | 'spongeNI' | 'lianaGenie3' | 'lianaPathway') {
    const input = evt.target as HTMLInputElement;
    const file = input.files && input.files[0] ? input.files[0] : undefined;
    if (!file) return;

    switch (type) {
      case 'spatial': this.spatialFile = file; break;
      case 'singleCell': this.singleCellFile = file; break;
      case 'multiome': this.multiomeFile = file; break;
      case 'genie3Net': this.genie3NetFile = file; break;
      case 'spongeNA': this.spongeNAFile = file; break;
      case 'spongeNI': this.spongeNIFile = file; break;
      case 'lianaGenie3': this.lianaGenie3File = file; break;
      case 'lianaPathway': this.lianaPathwayFile = file; break;
    }
  }

  // ---------- validation ----------
  private networkUploadsOk(): boolean {
    // VIPER needs Genie3
    const viper = !!this.form.get('network.algorithms.viper')?.value;
    if (viper && !this.genie3NetFile) return false;

    // AUCell/GSVA/ssGSEA need Genie3 OR both SPONGE files
    if (this.hasTriadSelected()) {
      const hasGenie3 = !!this.genie3NetFile;
      const hasSponge = !!this.spongeNAFile && !!this.spongeNIFile;
      if (!hasGenie3 && !hasSponge) return false;
    }
    return true;
  }

  canSubmit(): boolean {
    const spatialOk = !!this.spatialFile;

    const needsSingleCell = !!this.form.get('useTangram')!.value;
    const singleCellOk = !needsSingleCell || !!this.singleCellFile;

    const needsMultiome = !!this.form.get('useTangramMultiome')!.value;
    const multiomeOK = !needsMultiome || !!this.multiomeFile;

    const networkOn = !!this.form.get('scores.networkScores')?.value;
    const networkOk = !networkOn || this.networkUploadsOk();

    return spatialOk && singleCellOk &&  multiomeOK  && networkOk && !this.uploading;
  }

  // ---------- submit ----------
  async onSubmit() {
    if (!this.canSubmit()) return;

    this.errorMsg = '';
    this.resultJsonUrl = null;
    this.serverPayload = null;
    this.uploading = true;
    this.uploadProgress = 0;

    const fd = this.buildFormData();

    this.http.post<{ ok: boolean; json_url: string; json_filename: string; payload: any }>(
      '/api/upload',
      fd,
      { observe: 'events', reportProgress: true, withCredentials: true }
    ).subscribe({
      next: (evt: HttpEvent<any>) => {
        if (evt.type === HttpEventType.UploadProgress && evt.total) {
          this.uploadProgress = Math.round(100 * (evt.loaded / evt.total));
        }
        if (evt.type === HttpEventType.Response) {
          const body = evt.body || {};
          this.resultJsonUrl = body.json_url || null;
          this.serverPayload = body.payload || null;
          this.uploadProgress = 100;
          this.uploading = false;

          console.log('Upload finished, output_files', body.output_files)

          const geojsonPath = body.output_files?.geojsonPath; // e.g., "uploads/alice/results/hexagons.geojson"
          const parts = geojsonPath?.split('/');
          const user = parts?.at(-3);
          const subdir = parts?.at(-2);
          const filename = parts?.at(-1);

          console.log('filepath', geojsonPath)

          // --- update paths here ---
          const newPaths = {
            adataPath: body.output_files?.adataPath,
            genieFiltPath:  body.output_files?.genieFiltPath,
            spongeFiltPath: body.output_files?.spongeFiltPath,
            hexagonPath: `/api/hexagon/${user}/${subdir}/${filename}`, 
          };

          console.log('New Paths', newPaths)

          this.pathsService.updatePaths(newPaths);
        }
      },
      error: (err) => {
        this.errorMsg = err?.error?.detail || err.message || String(err);
        this.uploading = false;
      }
    });
  }

  private buildFormData(): FormData {
    const fd = new FormData();

    // --- core
    fd.append('email', this.form.value.email);
    fd.append('dataset', this.form.value.dataset);


    // --- spatial
    if (this.spatialFile) fd.append('spatial_h5ad', this.spatialFile);
    fd.append('spatial_normalization', String(this.form.value.spatialOptions.normalization));
    fd.append('spatial_filtering', String(this.form.value.spatialOptions.filtering));

    // --- tangram
    fd.append('use_tangram', String(this.form.value.useTangram));
    if (this.form.value.useTangram && this.singleCellFile) {
      fd.append('single_cell_h5ad', this.singleCellFile);
      fd.append('singlecell_filtering', String(this.form.value.tangram.filterSingleCell));
      fd.append('singlecell_normalization', String(this.form.value.tangram.normalizeSingleCell));
    }

     // --- multiome
    fd.append('use_multiome', String(this.form.value.useTangramMultiome));
    if (this.form.value.useTangramMultiome && this.multiomeFile) {
      fd.append('multiome_rds', this.multiomeFile);
    }

    // --- scores toggles
    const scores = this.form.value.scores;
    fd.append('score_network', String(scores.networkScores));
    fd.append('score_squidpy', String(scores.squidpy));
    fd.append('score_liana_plus', String(scores.lianaPlus));
    if (this.form.value.useTangramMultiome && this.multiomeFile) {
      fd.append('score_chromVar', String(scores.chromVAR));
      fd.append('score_diffMotifActivity', String(scores.diffMotifActivity));
      fd.append('score_motifEnrichment', String(scores.motifEnrichment));
      fd.append('score_FootprintingBias', String(scores.FootprintingBias));
      fd.append('genome', this.form.value.genome);
    }

    // --- network details
    if (scores.networkScores) {
      const alg = this.form.value.network.algorithms;
      fd.append('alg_viper', String(alg.viper));
      fd.append('alg_aucell', String(alg.aucell));
      fd.append('alg_gsva', String(alg.gsva));
      fd.append('alg_ssgsea', String(alg.ssgsea));

      // SPONGEeffects-Parameter nur wenn du sie zeigst (usingSponge==true)
      if (this.hasTriadSelected() && this.usingSponge) {
        const p = this.form.value.network.params;
        fd.append('net_m_score_threshold', p.mScoreThreshold ?? '');
        fd.append('net_p_adjust', p.pAdjust ?? '');
        fd.append('net_ensembl_id_col', p.ensemblId ?? '');
        fd.append('net_feature_col', p.featureCol ?? '');
        fd.append('net_rna_types', p.rnaTypes ?? '');
        fd.append('net_max_modules', p.maxModules ?? '');
      }

      // Genie3 parameter nur wenn Datei da
      if (this.genie3NetFile) {
        const g = this.form.value.network.genie3Params ?? {};
        fd.append('genie3_top_n_weights', String(g.topNWeights ?? 100000));
        fd.append('genie3_n_regulatory_genes', String(g.nRegulatoryGenes ?? 20));
        fd.append('genie3_n_regulons', String(g.nRegulons ?? 20));
        fd.append('genie3_network', this.genie3NetFile); // Datei selbst
      }

      // SPONGE Dateien (falls hochgeladen)
      if (this.spongeNAFile) fd.append('sponge_networkanalysis', this.spongeNAFile);
      if (this.spongeNIFile) fd.append('sponge_networkinteractions', this.spongeNIFile);
    }

    // --- squidpy (mehrere Methoden via Checkboxen)
    if (scores.squidpy) {
      const sq = this.form.value.squidpy;
      const m = sq.methods;
      if (m.moranI) {
        fd.append('squidpy_moranI', 'true');
        if (sq.moranI.nPerms !== null && sq.moranI.nPerms !== undefined && sq.moranI.nPerms !== '') {
          fd.append('squidpy_moranI_n_perms', String(sq.moranI.nPerms));
        }
        fd.append('squidpy_moranI_two_tailed', String(!!sq.moranI.twoTailed));
        fd.append('squidpy_moranI_corr_method', sq.moranI.corrMethod ?? '');
      }
      if (m.gearyC) {
        fd.append('squidpy_gearyC', 'true');
        if (sq.gearyC.nPerms !== null && sq.gearyC.nPerms !== undefined && sq.gearyC.nPerms !== '') {
          fd.append('squidpy_gearyC_n_perms', String(sq.gearyC.nPerms));
        }
        fd.append('squidpy_gearyC_two_tailed', String(!!sq.gearyC.twoTailed));
        fd.append('squidpy_gearyC_corr_method', sq.gearyC.corrMethod ?? '');
      }
      if (m.centrality_score) {
        fd.append('squidpy_centrality_score', 'true');
        fd.append('squidpy_centrality_score_cluster_key', sq.centrality_score.clusterKey ?? '');
      }
      if (m.co_occurrence) {
        fd.append('squidpy_co_occurrence', 'true');
        fd.append('squidpy_co_occurrence_cluster_key', sq.co_occurrence.clusterKey ?? '');
        fd.append('squidpy_co_occurrence_interval', sq.co_occurrence.interval ?? '');
        if (sq.co_occurrence.nSplits !== null && sq.co_occurrence.nSplits !== undefined && sq.co_occurrence.nSplits !== '') {
          fd.append('squidpy_co_occurrence_n_splits', String(sq.gearyC.nPerms));
        }
      }
      if (m.neighborhood_enrichment) {
        fd.append('squidpy_neighborhood_enrichment', 'true');
        fd.append('squidpy_neighborhood_enrichment_cluster_key', sq.neighborhood_enrichment.clusterKey ?? '');
        fd.append('squidpy_neighborhood_enrichment_library_key', sq.neighborhood_enrichment.libraryKey ?? '');
        fd.append('squidpy_neighborhood_enrichment_n_perms', sq.neighborhood_enrichment.nPerms ?? '');
      }
    }

    // --- LIANA+
    if (scores.lianaPlus) {
      const l = this.form.value.liana;
      fd.append('liana_composition_column', l.compositionColumn ?? '');
      if (this.lianaGenie3File)  fd.append('liana_genie3_network', this.lianaGenie3File);
      if (this.lianaPathwayFile) fd.append('liana_pathway_network', this.lianaPathwayFile);
    }

    return fd;
  }



  onReset() {
    this.form.reset({
      email: '',
      dataset: 'Visium',
      spatialOptions: { normalization: false, filtering: false },
      useTangram: false,
      tangram: { filterSingleCell: false, normalizeSingleCell: false },
      scores: { networkScores: false, squidpy: false, lianaPlus: false },
      network: {
        algorithms: { viper: false, aucell: false, gsva: false, ssgsea: false },
        params: {
          mScoreThreshold: '0.1',
          pAdjust: '0.05',
          ensemblId: 'ensembl_id',
          featureCol: 'feature_type',
          rnaTypes: 'lncRNA,protein_coding',
          maxModules: '20',
        },
        genie3Params: { topNWeights: 100000, nRegulatoryGenes: 20, nRegulons: 20 },
      },
      squidpy: {
        methods: {
          moranI: false,
          gearyC: false,
          centrality_score: false,
          co_occurrence: false,
          neighborhood_enrichment: false,
        },
        moranI: { nPerms: '', twoTailed: false, corrMethod: 'fdr_bh' },
        gearyC: { nPerms: '', twoTailed: false, corrMethod: 'fdr_bh' },
        centrality_score: { clusterKey: 'leiden' },
        co_occurrence: { clusterKey: 'leiden', interval: '50', nSplits: '' },
        neighborhood_enrichment: { clusterKey: 'leiden', libraryKey: '', nPerms: '1000' },
      },

      liana: { compositionColumn: 'tangram_ct_pred' },
    });

    // clear files
    this.spatialFile = undefined;
    this.singleCellFile = undefined;
    this.genie3NetFile = undefined;
    this.spongeNAFile = undefined;
    this.spongeNIFile = undefined;
    this.lianaGenie3File = undefined;
    this.lianaPathwayFile = undefined;


    this.errorMsg = '';
    this.uploadProgress = 0;
  }
}
