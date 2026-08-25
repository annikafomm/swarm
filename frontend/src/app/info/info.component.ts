import { Component, OnDestroy } from '@angular/core';
import { Location } from '@angular/common';
import { CommonModule } from '@angular/common';
import { HttpClient } from '@angular/common/http';
import { Subject } from 'rxjs';
import { takeUntil } from 'rxjs/operators';
import { DatasetService } from '../datasets.service';
import { Dataset } from '../datasets.service';
import { SessionService } from '../session.service';
import { GeoDataService } from '../geo-data.service';
import { MatButtonModule } from '@angular/material/button';
import { MatIconModule } from '@angular/material/icon';
import { MatDialog, MatDialogModule } from '@angular/material/dialog';
import { PrivacyDialogComponent } from '../privacy-dialog/privacy-dialog.component';

@Component({
  selector: 'app-info',
  standalone: true,
  imports: [CommonModule, MatButtonModule, MatIconModule, MatDialogModule],
  templateUrl: './info.component.html',
  styleUrl: './info.component.scss'
})
export class InfoComponent implements OnDestroy {
  selectedDataset: Dataset | null = null;
  /**
   * The GeoJSON currently loaded by the map, if any. Used only to report which optional
   * `meta` blocks this dataset actually ended up with — the catalog of possible keys is
   * static, but presence is per-dataset.
   */
  private geoMeta: Record<string, any> | null = null;
  private destroy$ = new Subject<void>();

  constructor(
    private location: Location,
    private datasetService: DatasetService,
    private http: HttpClient,
    private sessionService: SessionService,
    private geoDataService: GeoDataService,
    private dialog: MatDialog
  ) {
    // Subscribe to dataset selection changes
    this.datasetService.selectedDataset$
      .pipe(takeUntil(this.destroy$))
      .subscribe(dataset => {
        this.selectedDataset = dataset;
      });

    this.geoDataService.data$
      .pipe(takeUntil(this.destroy$))
      .subscribe(data => {
        this.geoMeta = data?.meta ?? null;
      });
  }

  goBack() {
    this.location.back();
  }

  openPrivacyPolicy(): void {
    this.dialog.open(PrivacyDialogComponent, {
      width: '850px',
      maxHeight: '90vh',
      autoFocus: 'first-button',
      panelClass: 'privacy-dialog-panel',
    });
  }

  /**
   * Download a file by path (works for binary formats like h5ad)
   * Uses the blob approach instead of opening in new tab
   */
  downloadFile(filePath: string): void {
    const encodedPath = encodeURIComponent(filePath);
    const downloadUrl = `${this.sessionService.apiUrl}/api/download/${encodedPath}`;

    // Fetch the file as a blob
    this.http.get(downloadUrl, {
      responseType: 'blob',
      withCredentials: true
    }).subscribe({
      next: (blob) => {
        // Create a temporary link element
        const link = document.createElement('a');
        const url = window.URL.createObjectURL(blob);
        link.href = url;

        // Extract filename from path
        const fileName = filePath.split('/').pop() || 'download';
        link.download = fileName;

        // Trigger download
        document.body.appendChild(link);
        link.click();

        // Cleanup
        document.body.removeChild(link);
        window.URL.revokeObjectURL(url);
      },
      error: (err) => {
        console.error(`Failed to download file: ${filePath}`, err);
        alert('Failed to download file');
      }
    });
  }

  /** Whether the full per-dataset inventory below the summary rows is expanded. */
  public detailsExpanded = false;
  /** When true, catalog entries this dataset doesn't have are filtered out of the inventory. */
  public hideUnavailable = false;

  toggleDetails(): void {
    this.detailsExpanded = !this.detailsExpanded;
  }

  toggleHideUnavailable(): void {
    this.hideUnavailable = !this.hideUnavailable;
  }

  /**
   * Every analysis stage the pipeline can record on a dataset, each reported as explicitly
   * enabled or disabled. This is a fixed catalog rather than a list of what happens to be
   * switched on, so a reader can tell "this dataset has no footprinting" apart from "this
   * build doesn't know about footprinting" — the two look identical in a present-only list.
   */
  getAnalyses(): Array<{ name: string; enabled: boolean; requirement: string; note: string }> {
    const ds = this.selectedDataset;
    return [
      {
        name: 'Tangram cell-type mapping', enabled: !!ds?.use_tangram, requirement: 'optional',
        note: 'Needs a single-cell reference. Produces per-spot cell-type composition.',
      },
      {
        name: 'Multiome (RNA + ATAC)', enabled: !!ds?.use_multiome, requirement: 'optional',
        note: 'Needs a 10x Multiome .rds. Unlocks all ATAC analyses below.',
      },
      {
        name: 'ChromVAR motif activity', enabled: !!ds?.use_chromvar, requirement: 'optional',
        note: 'Needs the Multiome .rds but NOT the ATAC fragments file. Spatially projected through the Tangram mapping. Also forces the two autocorrelations below.',
      },
      {
        name: 'Differential motif activity', enabled: !!ds?.use_differential_motif_activity, requirement: 'optional',
        note: 'Requires ChromVAR. Ranks motifs per group against the rest, using a spot-level grouping column.',
      },
      {
        name: 'Footprint analysis', enabled: !!ds?.use_footprinting, requirement: 'optional',
        note: 'Requires Multiome plus an ATAC fragments file and its .tbi index. Genome-wide footprints are spatial; the gene-specific footprints behind GRN evaluation are computed on dissociated clusters instead, not in space.',
      },
      {
        name: "Moran's I spatial autocorrelation", enabled: !!ds?.use_moranI, requirement: 'optional',
        note: 'On chromVAR motif scores. Near +1 clustered, 0 random, -1 dispersed. Implied by ChromVAR motif activity. Statistics stay motif-level even when shown per transcription factor, so one motif mapping to several TFs repeats its value.',
      },
      {
        name: "Geary's C spatial autocorrelation", enabled: !!ds?.use_gearyC, requirement: 'optional',
        note: 'On chromVAR motif scores. Runs 0-2 on the reverse scale to Moran\'s I: below 1 is positive autocorrelation, about 1 is random, above 1 is locally heterogeneous.',
      },
    ].filter((entry) => !this.hideUnavailable || entry.enabled);
  }

  /**
   * Catalog of the AnnData keys the pipeline reads or writes. Presence cannot be checked from
   * the browser (the .h5ad is never shipped to the client), so this is reference documentation
   * — it replaces the hand-maintained "Format Requirements" prose, which had drifted from the
   * code in three places (see the notes on the co-occurrence and neighbourhood-enrichment rows).
   */
  getAnnDataSpec(): Array<{ group: string; entries: Array<{ key: string; requirement: string; note: string }> }> {
    return [
      {
        group: 'Core matrix and annotations',
        entries: [
          { key: 'X', requirement: 'required', note: 'Expression counts, spots x genes.' },
          { key: 'obs', requirement: 'required', note: 'Per-spot metadata. Every column is surfaced as a GeoJSON property, so clean it before upload.' },
          { key: 'var', requirement: 'required', note: 'Per-gene metadata. Gene symbols only.' },
          { key: 'obsm["spatial"]', requirement: 'required', note: 'Spot coordinates. Without these no hexagon geometry can be built.' },
          { key: 'obs["in_tissue"]', requirement: 'optional', note: 'Visium tissue mask. Spots with in_tissue == 0 are dropped from the GeoJSON entirely, not just dimmed.' },
          { key: 'obs["leiden"]', requirement: 'conditional', note: 'Required for every cluster-level score. Computed automatically if a Squidpy method needs it and it is absent. Labels must cast to int — a non-integer label silently empties the whole cluster-annotation block.' },
          { key: 'obs["cell type"]', requirement: 'conditional', note: 'Required by Tangram in cluster mode (the default). Multiome forces cell-wise mode instead. A coarse annotation here can collapse every spot to one predicted type.' },
        ],
      },
      {
        group: 'Spatial autocorrelation',
        entries: [
          { key: 'uns["moranI"]', requirement: 'optional', note: 'Columns: I, pval_norm, var_norm, pval_norm_fdr_bh. Verified against current output. Also picks the gene the map opens on — the top row after sorting by I, with Moran\'s I taking priority over Geary\'s C.' },
          { key: 'uns["gearyC"]', requirement: 'optional', note: 'Columns: C, pval_norm, var_norm, pval_norm_fdr_bh. Verified against current output.' },
        ],
      },
      {
        group: 'Cluster-level scores',
        entries: [
          { key: 'uns["leiden_centrality_scores"]', requirement: 'optional', note: 'Columns: degree_centrality, average_clustering, closeness_centrality.' },
          { key: 'uns["leiden_co_occurrence"]', requirement: 'optional', note: 'Keys "occ" (n_clusters x n_clusters x n_intervals) and "interval" — singular, despite older docs saying "intervals".' },
          { key: 'uns["leiden_nhood_enrichment"]', requirement: 'optional', note: 'Keys "zscore" and "count" — singular, despite older docs saying "counts". Only zscore is read.' },
        ],
      },
      {
        group: 'Regulatory scores',
        entries: [
          { key: 'obsm["aucell_scores_genie3" | "_sponge"]', requirement: 'optional', note: 'AUCell gene-set activity per spot.' },
          { key: 'obsm["spongeffects_GSVA_scores_*"]', requirement: 'optional', note: 'GSVA module scores.' },
          { key: 'obsm["spongeffects_ssGSEA_scores_*"]', requirement: 'optional', note: 'ssGSEA module scores.' },
          { key: 'obsm["viper_scores_genie3"]', requirement: 'optional', note: 'VIPER regulon activity. GENIE3 only — there is no SPONGE variant.' },
          { key: 'uns["genie_genesets"] / uns["sponge_genesets"]', requirement: 'optional', note: 'Gene-set membership backing the network subgraph views.' },
        ],
      },
      {
        group: 'LIANA+ / cell communication',
        entries: [
          { key: 'obsm["ligand_receptor_cosine_similarity"]', requirement: 'optional', note: 'Per-spot ligand-receptor co-expression.' },
          { key: 'obsm["ligand_receptor_NMF_factors"]', requirement: 'optional', note: 'NMF factorisation of the LR signal.' },
          { key: 'obsm["cell_comp_tf_activity_cosine_similarity"]', requirement: 'optional', note: 'TF activity vs cell-type composition. Needs Tangram output.' },
          { key: 'obsm["tf_activity_score_ulm"]', requirement: 'optional', note: 'TF activity (decoupler ULM).' },
          { key: 'obsm["pathway_activity_score_mlm"]', requirement: 'optional', note: 'Pathway activity (decoupler MLM).' },
          { key: 'uns["liana_columns"]', requirement: 'conditional', note: 'Column names for the two similarity matrices above; without it they cannot be labelled.' },
        ],
      },
      {
        group: 'Multiome / ATAC',
        entries: [
          { key: 'obsm["chromvar_spot_scores"]', requirement: 'optional', note: 'chromVAR deviations projected onto spots.' },
          { key: 'uns["chromvar_motifs"]', requirement: 'conditional', note: 'Motif names for the matrix above; both are needed to reconstruct it.' },
          { key: 'uns["chromvar_moranI"] / ["chromvar_gearyC"]', requirement: 'optional', note: 'Spatial autocorrelation of motif activity.' },
          { key: 'uns["diff_motif_activity_top_motifs"]', requirement: 'optional', note: 'Per-comparison differential motif results.' },
          { key: 'uns[<grn_evaluation_name>]["peak_stats" | "motif_stats"]', requirement: 'optional', note: 'GRN evaluation tables, nested under the evaluation name. A flat uns["peak_stats"] is legacy and loses that name.' },
        ],
      },
    ].map(({ group, entries }) => ({ group, entries }));
  }

  /**
   * Catalog of the generated GeoJSON's shape. Unlike the AnnData spec, presence here is
   * checked live against the GeoJSON the map actually loaded for this dataset.
   */
  getGeoJsonSpec(): Array<{
    group: string;
    entries: Array<{ key: string; requirement: string; note: string; present: boolean | null }>;
  }> {
    const meta = this.geoMeta;
    const has = (key: string): boolean | null => {
      if (!meta) return null;
      const value = meta[key];
      if (value === undefined || value === null) return false;
      if (typeof value === 'object') return Object.keys(value).length > 0;
      return true;
    };

    const groups = [
      {
        group: 'Document structure',
        entries: [
          { key: 'type: "FeatureCollection"', requirement: 'required', note: 'One Feature per spot, each a hexagon Polygon.', present: meta ? true : null },
          { key: 'features[].geometry', requirement: 'required', note: 'Hexagon derived from obsm["spatial"] with the radius/scale options.', present: meta ? true : null },
          { key: 'features[].properties', requirement: 'required', note: 'Every obs column, plus the view-backing properties below.', present: meta ? true : null },
        ],
      },
      {
        group: 'Feature properties backing each view',
        entries: [
          { key: 'barcode', requirement: 'required', note: 'Spot identity; used as the D3 join key.', present: meta ? true : null },
          { key: 'leiden', requirement: 'conditional', note: 'Cluster id. Required by Cluster Information and Co-occurrence.', present: meta ? true : null },
          { key: 'regulatory_scores', requirement: 'optional', note: 'First available GENIE3/SPONGE score.', present: meta ? true : null },
          { key: 'chromvar_total_sum, motif_<group>', requirement: 'optional', note: 'chromVAR summaries per motif group.', present: meta ? true : null },
        ],
      },
      {
        group: 'meta blocks (checked live)',
        entries: [
          { key: 'data_type', requirement: 'required', note: 'visium or xenium; drives Xenium-specific rendering.', present: has('data_type') },
          { key: 'leiden_cluster_annotations', requirement: 'optional', note: 'Per-cluster centrality, co-occurrence and neighbourhood enrichment. Empty when Squidpy did not run.', present: has('leiden_cluster_annotations') },
          { key: 'interval', requirement: 'conditional', note: 'Co-occurrence distance-interval boundaries; needed to label the interval selector.', present: has('interval') },
          { key: 'moranI / gearyC', requirement: 'optional', note: 'Gene-level autocorrelation, backing Gene Expression.', present: has('moranI') || has('gearyC') },
          { key: 'grn_score_names', requirement: 'optional', note: 'Available regulatory score matrices.', present: has('grn_score_names') },
          { key: 'genie_genesets / sponge_genesets', requirement: 'optional', note: 'Gene-set dropdown contents.', present: has('genie_genesets') || has('sponge_genesets') },
          { key: 'ligand_receptor_global_scores', requirement: 'optional', note: 'Global LIANA+ ligand-receptor ranking.', present: has('ligand_receptor_global_scores') },
          { key: 'cell_comp_tf_activity_global_scores', requirement: 'optional', note: 'Global cell-composition TF activity ranking.', present: has('cell_comp_tf_activity_global_scores') },
          { key: 'tf_names / pathway_names / nmf_factors', requirement: 'optional', note: 'Column names populating the LIANA+ selectors.', present: has('tf_names') || has('pathway_names') || has('nmf_factors') },
          { key: 'chromvar_moranI / chromvar_gearyC', requirement: 'optional', note: 'Motif-level autocorrelation tables.', present: has('chromvar_moranI') || has('chromvar_gearyC') },
          { key: 'diff_motif_activity_top_motifs', requirement: 'optional', note: 'Differential motif activity per comparison.', present: has('diff_motif_activity_top_motifs') },
          { key: 'peak_stats / motif_stats', requirement: 'optional', note: 'GRN evaluation tables, keyed by evaluation name. Gates the GRN Evaluation tab.', present: has('peak_stats') || has('motif_stats') },
          { key: 'dgea', requirement: 'optional', note: 'Differential expression results per obs column.', present: has('dgea') },
        ],
      },
    ];

    return groups
      .map(({ group, entries }) => ({
        group,
        entries: entries.filter((entry) => !this.hideUnavailable || entry.present !== false),
      }))
      .filter(({ entries }) => entries.length > 0);
  }

  /**
   * The full catalog of files a dataset can carry, grouped by role. Entries keep their slot
   * with `path: undefined` when this dataset doesn't have them, so the inventory shows what
   * is missing as well as what is present.
   *
   * `geojson_path` is a `/api/geojson/<id>` endpoint rather than a file on disk, so it is
   * reported but never offered as a download — `/api/download` takes real paths only.
   */
  getFileGroups(): Array<{
    group: string;
    files: Array<{ name: string; key: string; path?: string; downloadable: boolean; requirement: string; note: string }>;
  }> {
    const ds = this.selectedDataset;
    const groups: Array<{
      group: string;
      files: Array<{ name: string; key: string; path?: string; requirement: string; note: string }>;
    }> = [
      {
        group: 'Core data',
        files: [
          {
            name: 'Primary AnnData', key: 'adata_path', path: ds?.adata_path, requirement: 'required',
            note: 'Not a file of its own — a pointer to whichever scored object this dataset is visualised from (see below).',
          },
          {
            name: 'Spatial scores AnnData', key: 'adata_st_scores_path', path: ds?.adata_st_scores_path, requirement: 'conditional',
            note: 'Scores computed directly on the spatial spots. Becomes the primary object when Tangram was not used.',
          },
          {
            name: 'Tangram-projected scores AnnData', key: 'adata_tg_scores_path', path: ds?.adata_tg_scores_path, requirement: 'conditional',
            note: 'Same scores computed on the Tangram-projected object. Becomes the primary object whenever Tangram or Multiome ran.',
          },
          {
            name: 'Tangram-projected AnnData (legacy alias)', key: 'tangram_adata_path', path: ds?.tangram_adata_path, requirement: 'optional',
            note: 'Documented as an alias for adata_tg_scores_path, so it points at projected output — not at the single-cell reference. The raw uploads are not retained as dataset fields at all.',
          },
          {
            name: 'Map geometry (GeoJSON)', key: 'geojson_path', path: ds?.geojson_path, requirement: 'required',
            note: 'Hexagon geometry plus all per-spot properties and the meta blocks. The map cannot render without it.',
          },
        ],
      },
      {
        group: 'Tangram mapping',
        files: [
          {
            name: 'Mapping matrix', key: 'adata_map_path', path: ds?.adata_map_path, requirement: 'conditional',
            note: 'Cells x spots probabilities. Required to project chromVAR scores onto spots.',
          },
          {
            name: 'Mapping matrix (CSV)', key: 'adata_map_X_csv_path', path: ds?.adata_map_X_csv_path, requirement: 'optional',
            note: 'Plain-text export of the matrix above.',
          },
          {
            name: 'Spot metadata (CSV)', key: 'adata_map_var_csv_path', path: ds?.adata_map_var_csv_path, requirement: 'optional',
            note: 'Spot-side annotations from the mapping step.',
          },
        ],
      },
      {
        group: 'Regulatory networks',
        files: [
          {
            name: 'GENIE3 network', key: 'genie_network_path', path: ds?.genie_network_path, requirement: 'conditional',
            note: 'Required for VIPER, and for AUCell/GSVA/ssGSEA unless SPONGE files are supplied instead.',
          },
          {
            name: 'SPONGE network', key: 'sponge_network_path', path: ds?.sponge_network_path, requirement: 'conditional',
            note: 'Alternative to GENIE3 for the AUCell/GSVA/ssGSEA trio.',
          },
          {
            name: 'TF regulatory graph', key: 'tf_graph_path', path: ds?.tf_graph_path, requirement: 'optional',
            note: 'Backs the GRN Evaluation network view.',
          },
        ],
      },
      {
        group: 'Multiome / ATAC',
        files: [
          {
            name: 'Seurat multiome object', key: 'global_motif_analysis_path', path: ds?.global_motif_analysis_path, requirement: 'optional',
            note: 'Full Seurat object carrying peaks, motifs and enrichment results.',
          },
          {
            name: 'Motif-to-TF mapping', key: 'motif_to_tf_csv_path', path: ds?.motif_to_tf_csv_path, requirement: 'optional',
            note: 'JASPAR motif id to transcription-factor name.',
          },
          {
            name: 'ChromVAR scores', key: 'chromvar_scores_csv_path', path: ds?.chromvar_scores_csv_path, requirement: 'optional',
            note: 'Per-cell motif deviations, before projection onto spots.',
          },
          {
            name: 'ChromVAR object (spot-level)', key: 'spot_obj_chromvar_path', path: ds?.spot_obj_chromvar_path, requirement: 'optional',
            note: 'Spot-projected chromVAR assay.',
          },
          {
            name: 'Footprints object (spot-level)', key: 'spot_obj_footprints_path', path: ds?.spot_obj_footprints_path, requirement: 'optional',
            note: 'Required for computing further footprints on demand.',
          },
          {
            name: 'Footprints object (dissociated)', key: 'dissociated_obj_footprints_path', path: ds?.dissociated_obj_footprints_path, requirement: 'optional',
            note: 'Cell-level counterpart of the object above.',
          },
        ],
      },
      {
        group: 'Pipeline',
        files: [
          {
            name: 'Pipeline log', key: 'calc_scores_log_path', path: ds?.calc_scores_log_path, requirement: 'optional',
            note: 'Timing and step-by-step record of the run that produced this dataset.',
          },
        ],
      },
    ];

    return groups
      .map(({ group, files }) => ({
        group,
        files: files
          .filter((file) => !this.hideUnavailable || !!file.path)
          .map((file) => ({
            ...file,
            downloadable: !!file.path && !file.path.startsWith('/api/'),
          })),
      }))
      .filter(({ files }) => files.length > 0);
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
