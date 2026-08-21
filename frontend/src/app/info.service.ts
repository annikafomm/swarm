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
    key: 'sample',
    title: 'Sample Information',
    icon: 'biotech',
    description:
      'Dataset-wide experimental metadata and biological origin details for this sample (e.g. sample ID, donor ID, treatment, condition, assay, organism, disease, tissue).',
  },
  {
    key: 'location',
    title: 'Location & Coordinates',
    icon: 'my_location',
    description:
      'Spatial position on the capture array grid (array row/column), pixel coordinates, and in-tissue status for this cell or spot.',
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

  sample: 'sample',
  sample_id: 'sample',
  donor: 'sample',
  donor_id: 'sample',
  patient: 'sample',
  patient_id: 'sample',
  treatment: 'sample',
  condition: 'sample',
  perturbation: 'sample',
  batch: 'sample',
  replicate: 'sample',
  is_primary_data: 'sample',
  suspension_type: 'sample',
  tissue_type: 'sample',

  assay: 'ontology',
  organism: 'ontology',
  tissue: 'ontology',
  disease: 'ontology',
  sex: 'ontology',
  development_stage: 'ontology',
  self_reported_ethnicity: 'ontology',
  assay_ontology_term_id: 'ontology',
  organism_ontology_term_id: 'ontology',
  tissue_ontology_term_id: 'ontology',
  cell_type_ontology_term_id: 'ontology',
  disease_ontology_term_id: 'ontology',
  development_stage_ontology_term_id: 'ontology',
  sex_ontology_term_id: 'ontology',
  self_reported_ethnicity_ontology_term_id: 'ontology',

  uniform_density: 'density',
  rna_count_based_density: 'density',

  array_row: 'location',
  array_col: 'location',
  spatial_x: 'location',
  spatial_y: 'location',
  x_original: 'location',
  y_original: 'location',
  _x: 'location',
  _y: 'location',
  x: 'location',
  y: 'location',
  in_tissue: 'location',
  observation_joinid: 'location',
};

const FIELD_INFO: Record<string, string> = {
  leiden: 'Cluster ID assigned by the Leiden community-detection algorithm on the cell/spot neighborhood graph. Cells in the same cluster share similar expression profiles.',
  leiden_centrality: 'Computed with Squidpy (sq.gr.centrality_scores): Spatial graph network centrality metrics (degree, closeness, and clustering coefficient) evaluating the spatial organization and connectivity of Leiden clusters.',
  degree_centrality: 'Computed with Squidpy (sq.gr.centrality_scores): Fraction of spatial neighborhood connections linked to cells in this cluster. Measures overall boundary surface and spatial graph connectivity.',
  average_clustering: 'Computed with Squidpy (sq.gr.centrality_scores): Average clustering coefficient of member cells, indicating how interconnected each cell\'s spatial neighbors are.',
  closeness_centrality: 'Computed with Squidpy (sq.gr.centrality_scores): Reciprocal of shortest spatial graph paths from this cluster to all other spots in the tissue. Higher values indicate a central spatial position.',
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
  spatial_x: 'Spatial X coordinate of this cell/spot in the tissue coordinate frame.',
  spatial_y: 'Spatial Y coordinate of this cell/spot in the tissue coordinate frame.',
  x_original: 'Original spatial X coordinate before registration or transformation.',
  y_original: 'Original spatial Y coordinate before registration or transformation.',
  _x: 'X coordinate or grid position of this cell/spot.',
  _y: 'Y coordinate or grid position of this cell/spot.',
  x: 'X coordinate of this cell/spot.',
  y: 'Y coordinate of this cell/spot.',
  in_tissue: 'Whether this spot was annotated as lying within the tissue section (1) or in background/off-tissue area (0).',
  sample: 'Identifier for the tissue sample this cell was collected from.',
  sample_id: 'Identifier for the tissue sample this cell was collected from.',
  donor: 'Identifier for the donor/individual this cell was collected from.',
  donor_id: 'Identifier for the donor/individual this cell was collected from.',

  barcode: 'Unique cell/spot barcode identifier from the sequencing library — the primary key for this cell.',
  centroid: "Pixel coordinates of this cell/spot's centroid in the original spatial image.",
};

const SCORE_INFO: Record<string, string> = {
  // Spatial Autocorrelation
  moranI:
    "Moran's I measures global spatial autocorrelation (range -1 to +1). Values near +1 indicate strong spatial clustering (neighboring spots have similar values); values near 0 indicate random spatial distribution; values near -1 indicate dispersion (checkerboard pattern).",
  gearyC:
    "Geary's C measures local spatial autocorrelation (range 0 to 2). Values < 1 indicate positive spatial clustering (similar neighbors); values = 1 indicate spatial randomness; values > 1 indicate negative autocorrelation (dissimilar neighbors).",
  pval_norm:
    'Analytical p-value derived under the null hypothesis of complete spatial randomness based on standard normal approximation.',
  pval_norm_fdr_bh:
    'False Discovery Rate (FDR) Benjamini-Hochberg adjusted p-value, accounting for multiple hypothesis testing across all evaluated features.',
  var_norm:
    'Theoretical variance of the spatial autocorrelation statistic under the null distribution of spatial randomness.',

  // Gene Regulatory Network Scores
  aucell:
    'AUCell evaluates whether a gene set or regulon is actively expressed in an individual spot by calculating the Area Under the Curve across the top-ranked expressed genes.',
  gsva:
    'SPONGeffects GSVA (Gene Set Variation Analysis): Computes non-parametric, unsupervised sample-wise gene set enrichment based on kernel estimation of cumulative distribution functions.',
  ssgsea:
    'SPONGeffects ssGSEA (single-sample GSEA): Calculates an enrichment score for each spot by integrating the empirical cumulative distribution functions of genes within versus outside the regulon.',
  viper:
    'VIPER (Virtual Inference of Protein-activity by Enriched Regulon analysis): Estimates protein/TF activity from the coordinated expression changes of its direct downstream target genes.',
  genie3_network:
    'GENIE3 regulatory network: Inferred using tree-based ensemble methods (Random Forest / ExtraTrees). Edge weights quantify the predictive importance of the regulator transcription factor for target gene expression.',
  sponge_network:
    'SPONGE ceRNA network: Inferred using multiple sensitivity correlation (mscor) to identify competing endogenous RNA interactions mediated by shared microRNA binding.',
  genie3_weight_cutoff:
    'Minimum edge weight threshold for GENIE3 interactions. Higher values filter for high-confidence, strongly predictive regulatory connections.',
  sponge_pvalue_cutoff:
    'Maximum adjusted p-value threshold for SPONGE ceRNA interactions. Lower values filter for statistically significant microRNA sponge pairs.',
  min_edges:
    'Controls the size of the displayed subgraph by selecting the top N most significant edges first, inferring the participating nodes, and reconstructing their interactions.',
  gprofiler:
    'Functional enrichment analysis in g:Profiler: Tests participating network genes against Gene Ontology (GO), KEGG, Reactome, and WikiPathways.',

  // LIANA+ Cellular Communication & Activity
  ligand_receptor:
    'LIANA+ bivariate spatial analysis measuring local intercellular communication between ligand- and receptor-expressing spots in the tissue microenvironment.',
  cosine_similarity:
    'Local Cosine Similarity: Quantifies directional alignment and co-expression intensity between ligand and receptor in neighboring spatial spots.',
  lr_p_value:
    'Spatial permutation p-value assessing whether observed local ligand-receptor co-localization is significantly greater than expected by random chance.',
  lr_category:
    'Spatial expression classification quadrant: High-High (both ligand and receptor enriched), High-Low, Low-High, or Low-Low in the local neighborhood.',
  nmf_factors:
    'Non-negative Matrix Factorization (NMF): Decomposes spatial cell-cell communication into low-dimensional latent spatial communication programs.',
  cell_comp_tf_activity:
    'LIANA+ bivariate similarity evaluating local spatial correlation between predicted cell-type abundance and transcription factor activity.',
  tf_activity_ulm:
    'decoupleR Univariate Linear Model (ULM) z-score estimating transcription factor activity from target gene expression weights in DoRothEA regulons.',
  pathway_activity_mlm:
    'decoupleR Multivariate Linear Model (MLM) score estimating signaling pathway activation from PROGENy responsive footprint genes.',

  // Leiden Clustering & Spatial Graph Metrics (Squidpy)
  degree_centrality:
    'Computed with Squidpy (sq.gr.centrality_scores): Fraction of all spatial neighborhood connections linked to spots in this cluster. Measures overall boundary surface and spatial graph connectivity.',
  average_clustering:
    'Computed with Squidpy (sq.gr.centrality_scores): Average clustering coefficient of member spots, indicating how interconnected each spot\'s spatial neighbors are in the neighborhood graph.',
  closeness_centrality:
    'Computed with Squidpy (sq.gr.centrality_scores): Reciprocal sum of shortest spatial graph paths from this cluster to all other spots in the tissue. Higher values indicate a central spatial position.',
  leiden_centrality:
    'Computed with Squidpy (sq.gr.centrality_scores): Graph-theoretic network centrality metrics (degree, closeness, and clustering coefficient) computed on the spatial neighborhood graph for each Leiden cluster.',
  nhood_enrichment:
    'Squidpy spatial permutation test assessing whether spots of cluster A neighbor spots of cluster B significantly more (or less) frequently than expected by random chance (z-score).',
  co_occurrence:
    'Conditional probability P(j|i) that a cell of cluster j occurs within a given Euclidean distance interval of a cell of cluster i, revealing multi-scale tissue organization.',

  // Multiome, ChromVAR & TF Footprinting
  chromvar_motif_activity:
    'Signac/chromVAR motif accessibility deviation z-score across JASPAR transcription factor binding sites, measuring chromatin accessibility enrichment relative to background peaks.',
  chromvar_moran:
    'Spatial autocorrelation of chromatin motif accessibility across the tissue using Moran\'s I.',
  diff_motif_activity:
    'Differential motif accessibility testing between spatial clusters (avg_diff effect size and FDR-adjusted p-value).',
  footprinting:
    'Tn5 transposase insertion bias-corrected chromatin accessibility profile around TF motif centers, showing cleavage protection at the protein binding site flanked by accessible cut sites.',
  peak_gene_links:
    'Statistical correlation between distal chromatin accessibility peaks and target gene transcription.',
  grn_filter_mode:
    'Filter mode for the tripartite TF-Peak-Gene network: Prior (strict, verified TF-motif-peak connections), Extended (intermediate confidence), or Full (all inferred links).',

  // Spatial Autocorrelation & Statistics
  getis_ord:
    'Getis-Ord General G measures global spatial clustering of high vs low values. A positive Z-score indicates spatial clustering of high values (hot spots), while a negative Z-score indicates spatial clustering of low values (cold spots).',
  moran_i:
    "Moran's I measures global spatial autocorrelation (range -1 to +1). Values near +1 indicate strong spatial clustering (neighboring spots have similar values); values near 0 indicate random spatial distribution; values near -1 indicate dispersion.",
  variance:
    'Sample variance quantifying the dispersion of values across all cells in this dataset.',
};

const COLUMN_INFO: Record<string, string> = {
  i: SCORE_INFO['moranI'],
  moran_i: SCORE_INFO['moran_i'],
  morani: SCORE_INFO['moran_i'],
  getis_ord: SCORE_INFO['getis_ord'],
  getisord: SCORE_INFO['getis_ord'],
  getis_ord_g: SCORE_INFO['getis_ord'],
  getis_ord_z: SCORE_INFO['getis_ord'],
  mean: 'Arithmetic mean across all valid spots.',
  variance: SCORE_INFO['variance'],
  var: SCORE_INFO['variance'],
  min: 'Minimum observed value across all spots.',
  max: 'Maximum observed value across all spots.',
  c: SCORE_INFO['gearyC'],
  pval_norm: SCORE_INFO['pval_norm'],
  pval_norm_fdr_bh: SCORE_INFO['pval_norm_fdr_bh'],
  var_norm: SCORE_INFO['var_norm'],
  logfc: SCORE_INFO['logFC'],
  p_val: 'Nominal p-value from the statistical significance test.',
  p_val_adj: 'Benjamini-Hochberg false discovery rate (FDR) adjusted p-value.',
  'p.adj': 'Benjamini-Hochberg false discovery rate (FDR) adjusted p-value.',
  'p.val': 'Nominal p-value from the statistical test.',
  weight: 'Predictive importance weight of the regulatory connection (GENIE3).',
  mscor: 'Multiple sensitivity correlation quantifying microRNA-mediated sponge interaction strength (SPONGE).',
  cor: 'Pearson correlation coefficient between expression profiles.',
  pcor: 'Partial correlation coefficient controlling for third-variable confounding.',
  degree: 'Number of direct interaction partners in the inferred regulatory network.',
  degree_centrality: 'Computed with Squidpy (sq.gr.centrality_scores): Fraction of spatial neighborhood connections linked to spots in this cluster.',
  average_clustering: 'Computed with Squidpy (sq.gr.centrality_scores): Average clustering coefficient of member spots in the spatial graph.',
  closeness_centrality: 'Computed with Squidpy (sq.gr.centrality_scores): Reciprocal sum of shortest spatial graph paths from this cluster to all other spots.',
  leiden_centrality: 'Computed with Squidpy (sq.gr.centrality_scores): Spatial graph network centrality metrics for Leiden clusters.',
  eigenvector: 'Eigenvector centrality measuring the influence of a node based on connections to other highly connected nodes.',
  betweenness: 'Betweenness centrality measuring how often a gene lies on the shortest path between other genes.',
  page_rank: 'PageRank score measuring network prominence and information flow.',
  cosine_similarity_mean: 'Mean local cosine similarity across all spatial spots for this interaction.',
  cosine_similarity_std: 'Standard deviation of local cosine similarity across spatial spots.',
  ligand_receptor_morans: "Moran's R spatial autocorrelation for this ligand-receptor interaction across the tissue.",
  avg_diff: 'Average difference in motif accessibility deviation score between the comparison groups.',
  pct_1: 'Percentage of cells in Group 1 with accessible signal for this feature.',
  pct_2: 'Percentage of cells in Group 2 with accessible signal for this feature.',
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
    if (/(?:^|_)(x|y)(?:_|$)/i.test(k)) {
      return 'location';
    }
    return 'other';
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

  /** Lookup detailed explanation for a score, metric, or analysis type. */
  getScoreInfo(key: string): string | null {
    if (!key) return null;
    const k = key.toLowerCase().replace(/[\s-]+/g, '_');
    return SCORE_INFO[k] || SCORE_INFO[key] || null;
  }

  /** Lookup detailed explanation for a table column header. */
  getColumnInfo(col: string): string | null {
    if (!col) return null;
    const clean = col.toLowerCase().trim();
    if (COLUMN_INFO[clean]) return COLUMN_INFO[clean];
    const underscored = clean.replace(/[\s.-]+/g, '_');
    if (COLUMN_INFO[underscored]) return COLUMN_INFO[underscored];
    return null;
  }
}
