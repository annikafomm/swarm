# SWARM Research Investigation & Paper Readiness Audit

## 1. Tool Capabilities Inventory

SWARM is an **interactive spatial multi-omics regulatory analysis web tool** that integrates:

| Layer | Method | What it computes | Output |
|:---|:---|:---|:---|
| **Transcriptional GRN** | GENIE3 → VIPER / AUCell | TF regulon activity per spatial spot | Per-spot regulatory scores, TF leaderboard |
| **Post-transcriptional ceRNA** | SPONGE → spongeffects (GSVA / ssGSEA / AUCell) | miRNA-mediated ceRNA module activity per spot | Per-spot sponge module scores |
| **Intercellular Signaling** | LIANA+ (bivariate local) | Ligand-receptor spatial communication, cell-composition × TF activity similarity | Per-spot LR scores, cosine similarity |
| **Pathway Activity** | decoupleR / PROGENy MLM | Signaling pathway activation per spot | Per-spot pathway scores |
| **Cell Deconvolution** | Tangram (+ cell2location) | Cell type abundances per spatial spot | 60+ cell type fractions (Heart) |
| **Spatial Statistics** | Squidpy | Moran's I, Geary's C, centrality, co-occurrence, nhood enrichment | Per-gene / per-cluster spatial metrics |
| **Epigenetic Validation** | Signac / chromVAR / Footprinting | Motif accessibility, differential motifs, TF binding footprints | Per-spot chromVAR scores, footprint PDFs |
| **GRN Benchmarking** | Peak-TF intersection | Precision-Recall of inferred GRN edges vs. ATAC peaks | Peak stats, motif stats, TF-Peak-Gene graph |
| **DGEA** | On-demand | Differential gene expression between arbitrary clusters | Volcano / ranked gene list |

### Unique Value Proposition
The tool's core novelty is **simultaneous interactive overlay** of all these layers on the same spatial coordinate system, enabling workflows that cross analytical boundaries (e.g., TF activity → chromatin validation → spatial co-localization → cell-cell communication) within a single interface.

---

## 2. Available Datasets

| # | Dataset | Tissue | Type | Spots | Cell Types | Scores Computed | Key Properties |
|:--|:---|:---|:---|:---:|:---:|:---|:---|
| 1 | **Builtin Visium** | Breast (TNBC) | Visium | 1,448 | 8 (stereoscope) | GENIE3 (AUCell/VIPER/GSVA/ssGSEA), LIANA+ (LR, TF-comp, pathway), Squidpy (Moran/Geary/centrality/co-occ/nhood), gene expression | `regulatory_scores`, `gene_expression`, `tf_activity`, `pathway_activity`, `ligand_receptor`, `cell_comp_tf_activity_similarity` |
| 2 | **Builtin Multiome** | Breast (TNBC) | Visium + ATAC | 1,448 | 8 | chromVAR (4 motif groups: GATA, RUNX, MYC, NFκB), motif spatial scores | `chromvar_total_sum`, `motif_GATA_like`, `motif_RUNX_like`, `motif_MYC_cluster`, `motif_NFkB_like` |
| 3 | **Uploaded Visium** | Visual Cortex (mouse) | Visium | 2,179 | — | GENIE3, Squidpy full suite, gene expression | `regulatory_scores`, `gene_expression`, `pathway_activity`, `ligand_receptor_relationships` |
| 4 | **Heart Cell Atlas** | Heart (human) | Multiome | 3,153 | 65 subtypes | Tangram deconvolution, chromVAR (4 groups), differential motif activity (6 comparisons), GRN evaluation, footprints (13 TFs) | 65 cell-type abundance columns (aCM1-5, vCM1-5, FB1-6, EC1-10, etc.), QC metrics, `annotation_final`, `broad_celltype` |

### Data Quality Assessment

| Criterion | Status | Notes |
|:---|:---:|:---|
| **Spot count** | ⚠️ Moderate | 1,448–3,153 spots. Sufficient for spatial statistics but small by Visium standards (typical: 2,000–5,000). |
| **Cell type diversity** | ✅ Good | Heart dataset has 65 fine-grained subtypes with expert annotation (`annotation_final`). Breast has 8 coarse types. |
| **Multi-modal coverage** | ✅ Strong | Heart has RNA + ATAC (true multiome). Breast has RNA only but with simulated ATAC overlay. |
| **GRN availability** | ⚠️ Partial | Breast Visium has GENIE3. Heart **does not** have GENIE3/SPONGE networks. Visual cortex has GENIE3. |
| **Spatial autocorrelation** | ⚠️ Partial | Computed for Breast Visium and Visual Cortex. **Not computed for Heart** (config shows `squidpy: false`). |
| **LIANA+ scores** | ⚠️ Partial | Computed for Breast Visium. **Not computed for Heart** (config shows `liana_plus: false`). |
| **Metadata completeness** | ⚠️ **Issue** | Builtin GeoJSONs have **empty `metadata: {}`** — meaning Moran's I tables, co-occurrence matrices, centrality scores, and DGEA are likely not being served to the frontend despite being computed in the h5ad. The `visium_to_geojson.py` conversion may need to be re-run. |

> [!WARNING]
> **Critical gap**: The Heart Cell Atlas dataset — the most biologically rich — has **no GRN (GENIE3/SPONGE), no LIANA+, and no Squidpy scores computed**. It only has Tangram deconvolution, chromVAR, differential motif activity, and GRN evaluation. This severely limits the cross-layer analysis workflows on the dataset that would produce the most interesting findings.

---

## 3. Workflows That Can Lead to Insights

### Workflow A: **Spatial GRN ↔ Chromatin Validation Loop** (⭐ Core Novelty)
1. Rank TFs by AUCell/VIPER activity (Tab: Regulatory Scores)
2. Identify spatially structured TFs via Moran's I
3. Validate: Does this TF's motif show accessible chromatin in the same spatial domain? (Tab: ChromVAR)
4. Confirm physical binding with Tn5 footprinting (Tab: Footprints)
5. Check GRN edge precision against ATAC peaks (Tab: GRN Evaluation)

**Status**: ⚠️ Only partially possible. Breast dataset has GRN but no ATAC. Heart has ATAC but no GRN. **Need a single dataset with both.**

### Workflow B: **Signaling → Transcription Cascade**
1. Identify spatially localized ligand-receptor interactions (Tab: LR Relationships)
2. Check which TFs are active in the receiving niche (Tab: TF Activity / Regulatory Scores)
3. Correlate: Does the spatial domain of high LR score overlap with high TF activity?

**Status**: ✅ Fully possible on Breast Visium dataset.

### Workflow C: **Cell Type Niche ↔ Regulatory Program**
1. Color map by deconvolved cell type abundance (e.g., aCM1)
2. Identify spatial co-localization with another cell type (Co-occurrence panel)
3. Check if spatially co-localized cell types show coordinated TF activity or LR signaling

**Status**: ⚠️ Partially possible. Heart has deconvolution but no LR/TF scores. Breast has LR/TF but only 8 coarse cell types.

### Workflow D: **Spatial Regulatory Discovery** (⭐ Highest Paper Impact)
1. Run DGEA between spatial clusters
2. Identify differentially expressed genes that are also targets of a spatially structured TF regulon
3. Check if the TF's motif is differentially accessible in chromVAR between these clusters
4. Validate with footprinting

**Status**: ❌ Not currently possible on any single dataset. Requires GRN + ATAC + LIANA + Squidpy on the same dataset.

---

## 4. Research Questions

### RQ1: Can SWARM identify spatially structured TF regulons whose activity is validated by chromatin accessibility?
- **Testable?** Only if GRN + ATAC exist on the same dataset (currently they don't).
- **Fix**: Recompute Heart dataset with GENIE3 network, or add ATAC data to Breast Visium.

### RQ2: Does intercellular ligand-receptor communication spatially co-localize with downstream TF activation?
- **Testable?** ✅ Yes, on Breast Visium (has LIANA+ LR + TF activity + regulatory scores).
- **Example investigation**: Do spots with high TGFB1-TGFBR2 ligand-receptor score also show high SMAD3 TF activity?

### RQ3: Do spatial regulatory network hubs (high Moran's I TFs) correspond to known tissue-specific master regulators?
- **Testable?** ✅ Yes, on Breast Visium (has Moran's I + GENIE3 regulons).
- **Example**: Rank TFs by spatial autocorrelation of their AUCell score → compare against known breast cancer master regulators (ESR1, FOXA1, GATA3, etc.).

### RQ4: Can SWARM reveal regulatory circuits that are invisible to either spatial expression analysis OR non-spatial GRN inference alone?
- **Testable?** ✅ Partially, on Breast Visium.
- **The key insight**: A TF might have high AUCell score (non-spatial regulon activity) but LOW Moran's I (spatially dispersed) → it's globally active but not tissue-structured. Conversely, a gene might have high Moran's I (spatially clustered) but NOT be part of any regulon → it's spatially structured but not transcriptionally regulated. **Only the intersection — high AUCell + high spatial autocorrelation — reveals the spatially organized regulatory program.**

### RQ5: Does the multi-scale co-occurrence of cell types predict specific intercellular signaling axes?
- **Testable?** ✅ Yes, on Breast Visium (co-occurrence + LIANA+ LR).

---

## 5. What The Tool Still Needs

### Critical for Paper

| Priority | Issue | Fix |
|:---:|:---|:---|
| 🔴 | **Heart dataset missing GRN + LIANA + Squidpy** | Re-run `calc_scores.py` on Heart with `-liana -moranI -gearyC -centrality_scores -co_occurrence -nhood_enrichment` and add a GENIE3 network. |
| 🔴 | **Builtin GeoJSON metadata is empty** | Re-run `visium_to_geojson.py` for all builtin datasets to populate Moran's I, co-occurrence, centrality, and DGEA tables in the frontend. |
| 🟡 | **No single dataset has all layers** (GRN + ATAC + LIANA + Squidpy + Tangram) | Create or identify one "showcase" dataset with complete analysis. Heart is the best candidate if GENIE3 + LIANA are added. |

### Nice-to-have for Paper

| Priority | Issue |
|:---:|:---|
| 🟢 | Cross-layer correlation view (scatter plot: TF AUCell score vs. chromVAR motif deviation) |
| 🟢 | Export: download spatial regulatory findings as a structured report |
| 🟢 | Batch comparison: compare regulatory programs across multiple tissue sections |

---

## 6. Findings That Require SWARM (Cannot Be Found Easily Otherwise)

### The "Intersection Insight" Pattern
Traditional analysis pipelines compute each layer independently:
- **Scanpy/Squidpy alone**: Finds spatially variable genes (Moran's I) but doesn't know which ones are TF targets.
- **GENIE3/SCENIC alone**: Infers TF→target networks but has no spatial information.
- **LIANA+ alone**: Finds LR interactions but doesn't know if TFs downstream are active.
- **chromVAR alone**: Finds accessible motifs but doesn't link to TF expression activity.

**SWARM uniquely enables**: "Show me TFs whose regulon is active (AUCell > threshold) **AND** whose targets are spatially autocorrelated (Moran's I > 0.3) **AND** whose motif is accessible in the same cluster (chromVAR deviation > 0) **AND** whose niche contains cells sending the right ligand."

This 4-layer intersection is the tool's killer feature and the basis for a compelling paper narrative.

---

## 7. Recommended Next Steps

1. **Immediate**: Re-run the Heart Cell Atlas dataset through the full pipeline with GENIE3 + LIANA + Squidpy enabled. This gives you a single dataset with all layers.
2. **Immediate**: Re-generate GeoJSON files with populated metadata (Moran's I tables, co-occurrence, centrality) for all builtin datasets.
3. **Investigation**: On the Breast Visium dataset, systematically cross-reference:
   - Top 20 TFs by AUCell spatial autocorrelation (Moran's I)
   - Their regulon targets' spatial autocorrelation
   - LIANA+ LR scores in the same spatial domain
   - Pathway activity in the same spatial domain
4. **Write-up**: Frame the paper around the "intersection insight" pattern, using one concrete biological finding (e.g., a specific TF circuit in breast cancer or heart tissue) as the central case study.
