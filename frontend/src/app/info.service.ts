import { Injectable } from '@angular/core';

export interface CellPropertyCategory {
  key: string;
  title: string;
  icon: string;
  description: string;
}

/**
 * Ordered top-to-bottom in the Cell Information tab. Order reflects what a biologist
 * typically wants to see first: what/where the cell is, then QC, then scores.
 */
const CATEGORIES: CellPropertyCategory[] = [
  {
    key: 'identity',
    title: 'Cluster & Cell Type',
    icon: 'category',
    description:
      'Cluster assignment and cell-type calls for this cell: the Leiden cluster ID, predicted or reference-annotated cell type, and any pathologist or computational-reference annotations available for this dataset.',
  },
  {
    key: 'composition',
    title: 'Cell-Type Composition',
    icon: 'pie_chart',
    description:
      "Estimated contribution of each reference cell type to this spot, from deconvolution of its (potentially mixed) transcriptome. Values are generally fractions that sum to roughly 1 across all listed cell types — useful for spots that capture more than one cell.",
  },
  {
    key: 'qc',
    title: 'Quality Control',
    icon: 'fact_check',
    description:
      'Standard sequencing-depth and library-complexity metrics for this cell or spot. Use these to judge data quality — e.g. a high mitochondrial fraction or very low gene/count numbers can indicate a stressed, dying, or low-quality cell.',
  },
  {
    key: 'regulatory',
    title: 'Regulatory & Motif Activity',
    icon: 'hub',
    description:
      "Gene-regulatory and chromatin-accessibility activity for this cell: AUCell scores summarizing regulon/gene-set activity from the GENIE3 and SPONGE gene regulatory networks, and chromVar deviation scores for transcription-factor motif accessibility.",
  },
  {
    key: 'liveView',
    title: 'Current View Value',
    icon: 'visibility',
    description:
      "The raw value(s) currently driving the map's color scale for this cell, snapshotted the last time the view was updated (e.g. the selected gene's expression, or the selected regulatory/pathway score).",
  },
  {
    key: 'ontology',
    title: 'Ontology Metadata',
    icon: 'schema',
    description:
      'Standardized ontology terms describing this cell\'s assay, tissue, organism, developmental stage, sex and disease context, following the CELLxGENE schema (Cell Ontology, Uberon, EFO, …) for cross-dataset comparability.',
  },
  {
    key: 'density',
    title: 'Spatial Density',
    icon: 'blur_on',
    description:
      "Local spatial density around this cell — how tightly packed neighboring cells/spots are, either by simple proximity (uniform) or weighted by each neighbor's RNA content.",
  },
  {
    key: 'location',
    title: 'Location & Sample',
    icon: 'my_location',
    description:
      'Spatial position on the capture array and sample/donor identifiers for this cell — useful for distinguishing cells from different biological replicates or tracing a spot back to its slide coordinates.',
  },
  {
    key: 'other',
    title: 'Additional Properties',
    icon: 'more_horiz',
    description: "Other properties present in this dataset that don't fall into the categories above.",
  },
];

const EXACT_FIELD_CATEGORY: Record<string, string> = {
  leiden: 'identity',
  cell_type: 'identity',
  annotation_jc: 'identity',
  annotation_final: 'identity',
  broad_celltype: 'identity',
  broad_celltype_score: 'identity',
  broad_celltype_top2: 'identity',
  broad_celltype_margin: 'identity',
  pathologist_annotation: 'identity',
  stereoscope_annotation: 'identity',

  n_genes_by_counts: 'qc',
  total_counts: 'qc',
  log1p_total_counts: 'qc',
  n_counts: 'qc',
  n_genes: 'qc',
  mt_frac: 'qc',
  pct_counts_mt: 'qc',
  pct_counts_in_top_50_genes: 'qc',
  pct_counts_in_top_100_genes: 'qc',
  pct_counts_in_top_200_genes: 'qc',
  pct_counts_in_top_500_genes: 'qc',

  chromvar_total_sum: 'regulatory',
  aucell_genie3: 'regulatory',
  aucell_sponge: 'regulatory',

  regulatory_scores: 'liveView',
  gene_expression: 'liveView',
  tf_activity: 'liveView',
  pathway_activity: 'liveView',
  ligand_receptor_relationships: 'liveView',
  cell_comp_tf_activity_similarity: 'liveView',

  assay: 'ontology',
  organism: 'ontology',
  tissue: 'ontology',
  disease: 'ontology',
  sex: 'ontology',
  development_stage: 'ontology',
  assay_ontology_term_id: 'ontology',
  organism_ontology_term_id: 'ontology',
  tissue_ontology_term_id: 'ontology',
  cell_type_ontology_term_id: 'ontology',
  disease_ontology_term_id: 'ontology',
  development_stage_ontology_term_id: 'ontology',
  self_reported_ethnicity_ontology_term_id: 'ontology',

  uniform_density: 'density',
  rna_count_based_density: 'density',

  array_row: 'location',
  array_col: 'location',
  in_tissue: 'location',
  sample: 'location',
  sample_id: 'location',
  donor: 'location',
  donor_id: 'location',
};

const FIELD_INFO: Record<string, string> = {
  leiden: 'Cluster ID assigned by the Leiden community-detection algorithm on the cell/spot neighborhood graph. Cells in the same cluster share similar expression profiles.',
  cell_type: "This cell's type label — either predicted computationally or taken from the dataset's existing annotation.",
  annotation_jc: 'A curated cell-type/cluster annotation supplied with this dataset.',
  annotation_final: 'The final, curated cell-type annotation for this cell after review.',
  broad_celltype: 'Coarse-grained (broad) cell-type call for this cell, grouping finer subtypes together.',
  broad_celltype_score: 'Confidence score for the broad cell-type call — higher means a more confident assignment.',
  broad_celltype_top2: 'The second most likely broad cell type for this cell, useful for spotting ambiguous calls.',
  broad_celltype_margin: 'Difference in score between the top and second cell-type call — a small margin flags an ambiguous assignment.',
  pathologist_annotation: 'Region or cell-type annotation manually assigned by a pathologist reviewing the tissue image.',
  stereoscope_annotation: 'Cell-type annotation derived from Stereoscope deconvolution of this spot.',

  n_genes_by_counts: 'Number of distinct genes detected with non-zero counts in this cell.',
  total_counts: 'Total number of UMI/read counts detected in this cell — a measure of sequencing depth.',
  log1p_total_counts: 'Natural log(1 + total counts) — a variance-stabilized version of sequencing depth, often used in QC plots.',
  n_counts: 'Total UMI/read count for this cell (legacy naming, equivalent to total counts).',
  n_genes: 'Number of genes detected in this cell (legacy naming, equivalent to n_genes_by_counts).',
  mt_frac: 'Fraction of total counts coming from mitochondrial genes. High values often indicate stressed or dying cells and are a standard QC filter.',
  pct_counts_mt: 'Percentage of total counts coming from mitochondrial genes. High values often indicate stressed or dying cells and are a standard QC filter.',
  pct_counts_in_top_50_genes: 'Percentage of this cell\'s total counts coming from just its top 50 most highly expressed genes. High values can indicate low library complexity.',
  pct_counts_in_top_100_genes: 'Percentage of this cell\'s total counts coming from its top 100 most highly expressed genes.',
  pct_counts_in_top_200_genes: 'Percentage of this cell\'s total counts coming from its top 200 most highly expressed genes.',
  pct_counts_in_top_500_genes: 'Percentage of this cell\'s total counts coming from its top 500 most highly expressed genes.',

  chromvar_total_sum: 'Summed chromVar deviation score across all transcription-factor motifs — a rough summary of overall regulatory/chromatin activity in this cell.',
  aucell_genie3: 'AUCell activity score for each GENIE3-inferred regulon/gene set in this cell: how strongly that gene set\'s expression is enriched relative to all genes.',
  aucell_sponge: 'AUCell activity score for each SPONGE-inferred regulon/gene set in this cell.',
  'motif_*': 'chromVar accessibility deviation score for this transcription-factor motif family — higher values indicate greater predicted regulatory activity of that motif in this cell.',

  regulatory_scores: 'Raw regulatory-network (GENIE3/SPONGE) AUCell score most recently used to color the map.',
  gene_expression: "Expression value of the currently selected gene, as last used to color the map.",
  tf_activity: 'Transcription-factor activity score (from motif/chromVar analysis) most recently used to color the map.',
  pathway_activity: 'Pathway activity score most recently used to color the map.',
  ligand_receptor_relationships: 'Ligand-receptor interaction score (LIANA+) most recently used to color the map.',
  cell_comp_tf_activity_similarity: "Similarity between this cell's composition and its transcription-factor activity profile, most recently used to color the map.",

  assay: 'Human-readable name of the experimental assay/technology used to generate this sample (e.g. 10x Visium).',
  organism: 'Organism this sample was taken from.',
  tissue: 'Tissue of origin for this sample.',
  disease: 'Disease/condition annotated for this sample, or "normal" if healthy.',
  sex: 'Sex of the donor this sample was taken from.',
  development_stage: 'Developmental stage of the donor at the time of sampling.',
  'assay_ontology_term_id': 'Experimental Factor Ontology (EFO) term ID for the assay — the machine-readable counterpart of "Assay".',
  'organism_ontology_term_id': 'NCBI Taxonomy ontology term ID for the organism — the machine-readable counterpart of "Organism".',
  'tissue_ontology_term_id': 'Uberon ontology term ID for the tissue — the machine-readable counterpart of "Tissue".',
  'cell_type_ontology_term_id': 'Cell Ontology (CL) term ID for the annotated cell type — the machine-readable counterpart of "Cell Type".',
  'disease_ontology_term_id': 'MONDO/PATO ontology term ID for the disease state — the machine-readable counterpart of "Disease".',
  'development_stage_ontology_term_id': 'Developmental-stage ontology term ID — the machine-readable counterpart of "Development Stage".',
  'self_reported_ethnicity_ontology_term_id': 'Ontology term ID for the donor\'s self-reported ethnicity, following the CELLxGENE schema.',
  '*_ontology_term_id': 'Machine-readable ontology term ID accompanying a human-readable metadata field, following the CELLxGENE schema.',

  uniform_density: 'Local cell density estimate using a uniform kernel over spatial neighbors, independent of expression.',
  rna_count_based_density: 'Local cell density estimate weighted by each neighbor\'s total RNA counts, so highly-expressing neighbors contribute more.',

  array_row: 'Row coordinate of this spot on the spatial capture array grid (e.g. a Visium slide), distinct from its physical pixel position.',
  array_col: 'Column coordinate of this spot on the spatial capture array grid (e.g. a Visium slide), distinct from its physical pixel position.',
  in_tissue: 'Whether this spot was annotated as lying within the tissue section (1) or in background/off-tissue area (0).',
  sample: 'Identifier for the tissue sample this cell was collected from.',
  sample_id: 'Identifier for the tissue sample this cell was collected from.',
  donor: 'Identifier for the donor/individual this cell was collected from.',
  donor_id: 'Identifier for the donor/individual this cell was collected from.',

  barcode: 'Unique cell/spot barcode identifier from the sequencing library — the primary key for this cell.',
  centroid: "Pixel coordinates of this cell/spot's centroid in the original spatial image.",
};

@Injectable({ providedIn: 'root' })
export class InfoService {
  readonly categories: CellPropertyCategory[] = CATEGORIES;

  /** Assigns a raw property key (and optionally its value) to a display category. */
  getCategoryForKey(key: string, value?: unknown): string {
    const k = key.toLowerCase();
    if (EXACT_FIELD_CATEGORY[k]) {
      return EXACT_FIELD_CATEGORY[k];
    }
    if (k.startsWith('motif_')) {
      return 'regulatory';
    }
    if (k.endsWith('_ontology_term_id')) {
      return 'ontology';
    }
    const isNumeric = typeof value === 'number' || (typeof value === 'string' && value.trim() !== '' && Number.isFinite(+value));
    return isNumeric ? 'composition' : 'other';
  }

  /** Short help text for a single field, or null if nothing specific is known about it. */
  getFieldInfo(key: string): string | null {
    const k = key.toLowerCase();
    if (FIELD_INFO[k]) {
      return FIELD_INFO[k];
    }
    if (k.startsWith('motif_')) {
      return FIELD_INFO['motif_*'];
    }
    if (k.endsWith('_ontology_term_id')) {
      return FIELD_INFO['*_ontology_term_id'];
    }
    return null;
  }
}
