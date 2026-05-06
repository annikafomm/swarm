
# TF GRN Evaluation Pipeline

## Overview

The **TF_GRN_evaluation** pipeline performs **Gene Regulatory Network (GRN) analysis** on multiome (ATAC + RNA) data. It identifies which transcription factors (TFs) regulate a given set of genes within specific cell clusters, by integrating:

1. **ATAC accessibility data** — which regulatory regions (peaks) are open in the chromatin
2. **RNA expression data** — which genes are expressed in each cell type
3. **Prior GRN knowledge** — predicted TF-gene regulatory relationships from databases
4. **Motif analysis** — identifying TF binding sites in accessible peaks

The pipeline is useful for understanding the regulatory logic driving gene expression in your cell type of interest, particularly for pathway-specific or disease-relevant gene modules.

## Three-Step Workflow

### Step 1: Get Gene Module (`01_get_gene_module.r`)
- **Input**: Gene Ontology (GO) ID or pathway identifier
- **Output**: `gene_module.csv` containing all genes in that GO term/pathway
- **Example**: GO term for "collagen fibril organization" → list of collagen genes + related genes

### Step 2: Get Prior GRN (`02_get_priorGRN.py`)
- **Input**: Gene module CSV from Step 1
- **Output**: `priorGRN.csv` — predicted TF-to-gene regulations from OmniPathDB database
- **Logic**: For each gene in your module, fetch known/predicted TFs that regulate it from literature/databases
- **Example**: For collagen genes, fetch all TFs known to regulate collagen production

### Step 3: Full GRN Evaluation (`03_GRN_evaluation.r`)
- **Input**: Multiome Seurat object, gene module, prior GRN, target cell cluster
- **Process**:

  #### 1. **Peak-Gene Linking** (`linking_peaks_for_module_cluster()`)
  Identifies which ATAC peaks are functionally linked to seed genes via a multi-step approach:
  - **GC normalization**: Computes GC content for each peak (for later background matching)
  - **Correlation-based linking**: Uses Signac's `LinkPeaks()` to find peaks with correlated accessibility & gene expression
  - **Distance filtering**: Keeps links within ±150 kb of gene TSS (removes very distant distal links)
  - **Proximal/Distal annotation**: Classifies each link as:
    - **Proximal**: 100–2000 bp upstream of TSS (promoter-bound)
    - **Distal**: Everything else (enhancer-bound)
  - **Cluster-specific filtering**: If analyzing a specific cluster, only considers cells in that cluster for correlation
  - **Output**: `res_links$object` with updated `Links()` metadata; `res_links$plot_links` histogram of link distances

  #### 2. **Peak Statistics** (`do_peak_stats()`)
  Computes per-peak, per-gene, per-cluster statistics to assess whether peaks are "active regulators":
  - **Accessibility test**: One-sided t-test comparing peak accessibility in target cluster vs. all other cells
    - Tests `H₀: accessibility(cluster) ≤ accessibility(other)` 
    - Outputs: `t-stat_highly.acc`, `p.value-t.test_highly.acc_BH` (BH-adjusted p-value)
  - **Peak clustering**: Groups peaks by GC content and length to create "comparable" background peaks for enrichment testing
  - **Zero-expression statistics**: For each (peak, gene, cluster) combination, computes:
    - `acc_not_0`: Fraction of cluster cells where peak is accessible
    - `expr_not_0`: Fraction of cluster cells where gene is expressed
    - `expr_not_0.given_acc`: P(gene expressed | peak accessible) in cluster
    - `FC.expr_given_acc`: Fold-enrichment of expression among accessible cells
    - Background versions (`.bg`, `.all`) for comparison across clusters
  - **Peak counts**: Tallies proximal vs. distal linked peaks per (gene, cluster)
  - **Filtering**: High-confidence peaks must pass thresholds on accessibility, expression correlation, and signal strength
  - **Output**: `object@misc$peak_stats` table; `object@misc$df.stats.filtered` high-confidence subset

  #### 3. **Motif Enrichment Analysis** (`do_local_motif_stats()`)
  Tests whether transcription factor binding motifs are enriched in peaks linked to seed genes:
  - **Motif scanning**: Uses JASPAR2024 & TFBSTools to identify all known TF motifs in linked peaks
  - **Proximal motif enrichment**:
    - Counts occurrences of each motif in proximal linked peaks
    - Compares to matched background: for each observed peak, randomly samples a GC-matched background peak and re-counts
    - t-test: `t.test(bg_motif_counts, mu = observed_count, alternative = "less")`
    - Outputs: `proximal.motif_count`, `proximal.background_count`, `log2FC.proximal`, `p_adjust.proximal`
  - **Distal motif enrichment**: Same logic applied to distal peaks
  - **Promoter motif scanning**: Scans 2000 bp upstream of gene TSS (±extended if proximal peaks exist) for motif binding sites
  - **Filtering**: Extends motifs for TFs in priorGRN by de nove motifs with significant enrichment (adjustable p-value & log2FC thresholds)
  - **Output**: `object@misc$motif_stats` table with all motif-level statistics

  #### 4. **Footprinting** (`add_motif_stats()` → `footprint_stats_test()`)
  Detects evidence of transcription factor binding via Tn5 insertion patterns in ATAC-seq:
  - **Principle**: Bound TFs protect DNA from Tn5 cleavage, creating a "footprint" (depleted insertions) around the binding site
  - **Foreground**: Tn5 insertion profile around motif sites in linked peaks (bias-corrected)
  - **Background**: Insertion profile around the same motifs in random genomic regions
  - **Score computation**: `footprint_score = mean(core insertions) - mean(flank insertions)`
    - More negative = stronger footprint (better evidence of TF binding)
    - Values near 0 = no footprint (inactive motif or no TF binding)
  - **Statistical test**: `t.test(background_scores, mu = foreground_score, alternative = "greater")`
    - Small p-value = foreground significantly more depleted than background = stronger footprint
  - **Output**: `footprint_score`, `footprint.p_value_adj`, `bg_footprint_mean`, coverage metrics (`left_flank_nonzero_positions`, `right_flank_nonzero_positions`)

  #### 5. **Network Construction** (`draw_subNetwork_for_context()`)
  Builds a directed gene regulatory network from all TF-gene-motif evidence:
  - **TF-motif mapping**: Maps detected motifs to known transcription factors via `add_motif_TF_mapping()`
    - Uses JASPAR2024 metadata to link each motif to its TF(s)
  - **Prior GRN overlap**: For each TF-gene pair:
    - Checks if relationship exists in prior GRN (from Step 2)
    - Colors green if TF is in priorGRN (supporting external evidence)
    - Colors blue if TF is not in priorGRN (novel prediction)
    - Colors red/grey if no linked peaks (conflicting or absent)
  - **Network graph construction**:
    - **Nodes**: All seed genes (target genes) + all TFs with footprinted motifs
    - **Edges**: TF → Gene with attributes:
      - `regulatorType`: proximal (solid line) vs. distal (dashed line)
      - `color`: indicates motif enrichment & footprinting strength
      - `priorTF`: boolean indicating prior GRN membership
  - **Scoring & filtering**: Edges ranked by combined evidence:
    - Motif enrichment (log2FC)
    - Footprint score (TF binding evidence)
    - Prior GRN support
  - **Output**: D3.js-compatible `graph.json` (nodes list + edges list with coordinates for visualization)
- **Outputs**:
  - `peak_stats.csv` — detailed peak-gene link statistics
  - `motif_stats.csv` — TF motif enrichment and footprint scores
  - `graph.json` — D3.js-compatible network visualization (TFs → genes)
  - `GRN_object.rds` — full Seurat object with all analysis results

## Running the Pipeline

### Quick Start (R)
```r
# Step 1: Get gene module from GO term
Rscript 01_get_gene_module.r "GO:0030414" "output/"
# Creates: output/gene_module.csv

# Step 2: Get prior GRN
python 02_get_priorGRN.py output/gene_module.csv output/priorGRN.csv

# Step 3: Run full GRN evaluation (set parameters in file or via command line)
Rscript 03_GRN_evaluation.r \
  --analysis_name "collagen_fibral_orga" \
  --gene_set_path "output/gene_module.csv" \
  --priorGRN_path "output/priorGRN.csv" \
  --global_analysis_out_path "global_motif_analysis.rds" \
  --cluster "Fibroblast" \
  --ident "cell_type" \
  --jobdir_dir "/jobs"
```

### Output Structure
```
/jobs/multiome/GRN_evaluation/{initial_upload|on_demand}/{analysis_name}/
├── peak_plots/
│   ├── umap.png                    # UMAP colored by cell type
│   ├── module_enrichment.pdf       # Module score enrichment test
│   └── links_distr.png             # Distribution of peak-gene links
├── peak_stats.csv                  # All peak-gene link statistics
├── motif_stats.csv                 # TF motif enrichment & footprinting
├── graph.json                      # D3.js network visualization
└── GRN_object.rds                  # Seurat object with detailed results
```

## Key Statistics Explained

See sections below for detailed definitions of all output columns.

---

# explainations of stats values

## explainations peak-stats 
#### Identity columns 
- 'cluster'
- 'peak'
- 'gene'
- 'annotation': seed gene, or module gene (from prior net), or marker gene
- 'regulatorType': proximal: 100-2000 upstream TSS, else distal

#### LinkPeaks stats: across all clusters
- 'signac.scores': The raw Signac link score from Links(object)$score
- 'signac.zscores': The Signac link z-score from Links(object)$zscore
- 'signac.pvalues': The Signac link p-value from Links(object)$pvalue

#### Cluster-specific peak accessibility test columns: conduct_stat_test() using test='t-test', test_activation=TRUE (is this peak more accessible in the target cluster than outside it?)
foreground: cells of cluster,  background: all other cells, alternative: greater
- t-stat_highly.acc: The one-sided t-test statistic for target-cluster accessibility > background accessibility
- p.value-t.test_highly.acc: The raw p-value from that one-sided t-test
- p.value-t.test_highly.acc_Bonf: Bonferroni-adjusted p-value across all unique (cluster, peak) tests
- p.value-t.test_highly.acc_BH: Benjamini–Hochberg adjusted p-value across all unique (cluster, peak) tests

#### Probability / zero-expression columns: zero_expression_stats() and calc_prob_stats(gene_expr, peak_acc)
(n_cells = number of cells in cluster
n_expr = # cells with gene_expr != 0
n_acc = # cells with peak_acc != 0
n_joint = # cells with gene_expr != 0 AND peak_acc != 0)

- expr_not_0 = n_expr / n_cells: Fraction of cells in the target cluster where the gene is expressed
- acc_not_0 = n_acc / n_cells: Fraction of cells in the target cluster where this peak is accessible
- expr_not_0.given_acc = n_joint / n_acc: Conditional probability that the gene is expressed among cells where this peak is accessible, within the target cluster: P(expr != 0 | acc != 0) in the cluster
- expr_and_acc_not_0 = n_joint / n_cells: Joint fraction of target-cluster cells with both gene expression and peak accessibility nonzero: P(expr != 0 AND acc != 0) in the cluster.
- FC.expr_given_acc = expr_not_0.given_acc / expr_not_0: Fold-enrichment of gene expression among accessible cells relative to the cluster baseline expression rate: So values above 1 mean expression is enriched among cells where the peak is open

#### .bg: Background: same peak, other clusters: 
These are computed on the same peak and same gene, but using all cells outside the target cluster.
.all: same peak, All cells:
These are the same metrics as above but computed using all cells for the same gene and same peak.

#### Per-seed peak counts
- 'promotersLinkedToSeed': Number of linked peaks for this (gene, cluster) whose regulatorType == "proximal"
- 'distalPeaksLinkedToSeed': Number of linked peaks for this (gene, cluster) whose regulatorType == "distal"

#### Convenience / derived filter columns
- 'acc_cells_cluster': Estimated number of accessible cells in the target cluster for this peak: acc_cells_cluster = acc_not_0 * cluster_size where cluster_size = table(Idents(object))[cluster]
- 'delta_expr_given_acc.same_peak_bg': Difference in conditional expression between target cluster and background, for the same peak: expr_not_0.given_acc - expr_not_0.given_acc.bg
(- 'delta_expr_given_acc.other_peaks.same_cluster': Difference between the target peak and the “other linked peaks of same gene” summary, inside the target cluster:
expr_not_0.given_acc - expr_not_0.given_acc.bg_other_peaks.same_cluster
'delta_expr_given_acc.other_peaks.all': Difference between the target peak’s across-all-cells conditional expression and the “other linked peaks of same gene” across-all-cells summary: expr_not_0.given_acc.all - expr_not_0.given_acc.bg_other_peaks.all)

#### Final filter flags
- 'pass_cluster_specific'
    - promoter requirement, if enabled: promotersLinkedToSeed > 0
    - p.value-t.test_highly.acc_BH < th
    - t-stat_highly.acc > cluster_t_min
    - acc_cells_cluster >= min.cells
    - expr_and_acc_not_0 >= cluster_expr_given_acc_min
    - expr_not_0.given_acc > expr_given_acc_th
    - FC.expr_given_acc > cluster_fc_min
    - delta_expr_given_acc.same_peak_bg >= cluster_delta_same_peak_bg_min
- 'pass_global'
    - promoter requirement, if enabled
    - signac.zscores >= global_signac_z_min
    - signac.pvalues < global_signac_p_cutoff
    - expr_and_acc_not_0.all >= global_expr_given_acc_min
    - expr_not_0.given_acc.all > expr_given_acc_th
    - FC.expr_given_acc.all > global_fc_min
- 'pass_any': TRUE if pass_cluster_specific | pass_global
- 'pass_type'


## motif stats: Motif enrichment and footprint scoring
#### Identity Columns
- cluster
- gene
- motif
- TF: Transcription factor that binds motif (if multiple, new row for each)

#### Motif enrichment in linked peaks: motif_enrichment_per_gene() / calculate_enrichments()
logic: 
1. find peaks linked to the gene in the chosen cluster
2. split them into proximal vs distal using regulatorType
3. count how often each motif occurs in those linked peaks
4. compare that to a random background made by replacing each linked peak with a sampled “comparable”(GC content and length) peak from the same peak meta-feature cluster. t-test for enrichment
-> for this: make clusters of comparable peaks, store in object$peaks@meta.features$cluster
##### Proximal motif columns
- proximal.motif_count: The number of occurrences of this motif in proximal peaks linked to the gene in this cluster
- proximal.background_count: The mean motif count in matched random background peak sets for the proximal linked peaks. For each background draw, every linked proximal peak is replaced by a sampled peak from the same object$peaks@meta.features$cluster, then motif counts are computed, and the mean over draws is stored here
- log2FC.proximal: log2(proximal.motif_count / proximal.background_count): log2 enrichment of proximal motif count over matched background. Positive values mean enrichment; zero means no enrichment; negative values mean depletion
- t_stat.proximal: The t-statistic from testing the proximal background distribution against the foreground proximal count using
t.test(col_bg, mu = value_fg, alternative = "less") in effect. Interpreted practically, small p-values support that the background mean is lower than the foreground count, i.e. the motif is enriched in the linked proximal peaks.
- p_value.proximal: The raw p-value for that proximal enrichment test.
- p_adjust.proximal: A simple multiplicity-adjusted proximal p-value, computed as raw p_value.proximal * number_of_motifs, then capped at 1. This is Bonferroni-like, not BH/FDR.

##### Distal motif columns: analog to proximal, but using distal peaks
- distal.motif_count
- distal.background_count
- log2FC.distal = log2(distal.motif_count / distal.background_count)
- t_stat.distal
- p_value.distal
- p_adjust.distal

##### Promoter motif column
- promoter.motif_count: The number of occurrences of this motif in the gene’s promoter sequence, regardless of whether there is an accessible linked peak there. This is computed by get_tf_bindingsites_in_region(), which scans the promoter DNA sequence with the motif PWM. The promoter length: 2000 bp upstream TSS, but it may be extended if the gene has proximal linked peaks farther upstream than 2 kb.

#### Footprint columns: add_motif_stats() -> footprint_stats_test() for each  (gene, cluster, motif)
idea: The foreground footprint is built from motif sites found in peaks linked to the gene; the background footprint distribution is built by repeatedly sampling the same number of motif sites of that motif from elsewhere and rescoring them
core = mean.footprint - mean.flanks
- difference = observed.insertions_normalized - expected.insertions 
- mean.flanks = mean(difference in left + right flank windows)
- mean.footprint = mean(difference in motif-core window)
more negative footprint_score means the motif center is more depleted relative to the flanks, which is the classic “footprint” pattern
values near 0 mean little difference between core and flanks
positive values would mean the core is more accessible than the flanks

observed.insertions_normalized: empirical insertion profile around the motif sites, normalized by its global mean (it is the mean observed Tn5 insertion signal at each relative position, normalized so the average across positions is 1)
expected.insertions: sequence-bias-based expected insertion profile, normalized by flank expectation  (comes from GetExpectedInsertion(), which extracts the DNA sequence around the motif regions, gets the assay’s Tn5 bias vector, and calls FindExpectedInsertions(). That helper computes the expected insertion profile from local sequence composition and Tn5 hexamer bias, then normalizes that expected profile by the mean of the flank positions.)

##### Footprint statistic columns
- footprint_scores: as explained above
- bg_size: how many bg scores where computed
- bg_footprint_mean: The mean of the background footprint scores across those sampled background footprints
- footprint.t_stat: t.test(x = background_scores, mu = foreground_score, alternative = "greater"): So the test is asking whether the background mean score is greater than the foreground score. Because more negative scores are stronger footprints in your setup, a small p-value supports the foreground having a stronger depletion footprint than background
- footprint.p_value: The raw p-value for that footprint-vs-background test.
- footprint.p_value_adj: A per-(gene, cluster) multiplicity-adjusted version of footprint.p_value. In add_motif_stats(), after all motif rows are combined, the code multiplies each raw footprint p-value by the number of motif rows for that same gene and cluster, then caps at 1. So again this is Bonferroni-like, not BH/FDR.

##### Footprint quality / coverage columns
- sd.flanks: The standard deviation of the flank difference values for the foreground footprint. This is a variability measure for the left/right flank regions around the motif
- bg_sd_mean: The mean of the flank standard deviations across the sampled background footprints
- left_flank_nonzero_positions: How many positions in the left flank had nonzero observed insertion counts in the foreground footprint pileup. It is computed from observed.insertions.counts_sum in the footprint plot data
- right_flank_nonzero_positions: same for right flank