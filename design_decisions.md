# SWARM Design Decisions & UX Options

## Overall structure

### Option 1: Big table 
Idea: aggregate according to features: 

- most tables have gene/gene product as row. We could aggregate them in a huge table --> clean 
- remaining extra tables: 
  - ligand receptor relationships: or we put them as a dropdown or sth into the big table
  - pathway activity: clearly seperate
  - DGEA: functional 
  - ChromVar spatial correlation? can't check rn 
  - Differential motif activity? can't check rn 
  - Footprints ? can't check rn 
  - GRN Evaluation & GRN evaluation on demand ? can't check rn 

### Another Layer: Global, local, differential: 
Most scores would be interesting on all of those layers: 
- Globally on the whole tissue area (this is what we mostly have currently)
- Locally subsetted to one cluster (leiden / celltype / etc): e.g. rank by absolute gene expression only whithin this cluster
- Differential: like the DGEA tab but also compare other scores than gene expression between clusters
- Spatial autocorrelation scores can also be computed Globally and Locally: locally within once cluster is e.g. interesting with gradients along a neural layer etc. I see it as a metric next to max, mean, var

### Option 2: Regulatory-First Architecture

**Problem**: The tool currently has 14 flat, unstructured tabs. Half of them are disabled or empty on datasets lacking ATAC or SPONGE, and generic features (QC, cell info, single-gene expression) overshadow SWARM's core novelty.

**Guiding Principle**: SWARM's primary scientific novelty is **mapping, scoring, and validating inferred gene regulatory networks (GENIE3), ceRNA crosstalk (SPONGE), and chromatin states (ATAC) on spatial coordinates**. Standard spatial tasks (cell typing, QC, single-gene inspection) serve as **supporting context**, not the hero experience.

---

#### The 4 Regulatory Pillar Tabs + Context Drawer

Instead of 14 flat tool tabs, the interface is organized into **4 core regulatory question tabs** with a persistent **Context & Inspector** drawer:

```
┌────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│  [ 🌲 1. Transcriptional GRNs ]  [ 🧽 2. ceRNA Networks ]  [ 🔄 3. Niche & Signaling ]  [ 🔓 4. Epigenetic Validation ]  │
└────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

##### 🌲 Tab 1: Transcriptional GRNs (GENIE3 & Regulons)
* **Goal**: Explore which transcription factors drive spatial gene programs and discover their downstream targets.
* **Contents**:
  - **Regulon Leaderboard**: Ranked TFs $\to$ 1-click project AUCell / ULM regulon score onto the spatial map.
  - **Target Gene Explorer**: Edge-weighted target gene list ($w_{ij}$) for the active TF.
  - **Cell Composition $\times$ TF Activity**: Sender cell type abundance $\to$ receiver TF activation.
* **Selection & Controls**: Score metric (`AUCell`, `Moran's I`, `Geary's C`, `ULM`), TF search, min regulon size filter, color scale palette & percentile clipping.

##### 🧽 Tab 2: Post-Transcriptional ceRNA Networks (SPONGE)
* **Goal**: Explore microRNA-mediated competitive endogenous RNA crosstalk and post-transcriptional spatial modules.
* **Contents**:
  - **Spatial ceRNA Module Activity**: spongeffects scores (`AUCell`, `GSVA`, `ssGSEA`) on spatial spots.
  - **ceRNA Interaction Matrix**: Significant gene-gene pairs ($mscor$, partial correlation, $p$-adj).
  - **Hub ceRNAs**: Most central sponging RNAs across the tissue.
* **Selection & Controls**: Scoring algorithm, $mscor$ cutoff slider, FDR $p$-value threshold.

##### 🔄 Tab 3: Spatial Signaling & Niches (LIANA+)
* **Goal**: Identify intercellular communication that triggers intracellular regulatory programs.
* **Contents**:
  - **Spatial Ligand-Receptor Table**: Ranked by Cosine Similarity, $p$-value, categories (e.g. *High L - High R*), global LR Moran's I.
  - **Spatial Co-occurrence & Neighborhoods (Squidpy)**: Co-occurrence probability curves over distance radius $r$, neighborhood enrichment heatmap.
* **Selection & Controls**: Source $\to$ Target cell type filters, category filter, distance scale slider ($r$).

##### 🔓 Tab 4: Epigenetic & Benchmark Validation (Multiome ATAC)
* **Goal**: Validate whether inferred GRNs reflect real, physical DNA-binding and accessible chromatin.
* **Contents**:
  - **chromVAR Motif Accessibility**: Spatial deviation z-scores (JASPAR 2024 CORE) + Moran's I / Geary's C.
  - **Differential Motif Activity**: Cluster-specific enriched motifs ($\Delta\text{accessibility}$, FDR).
  - **Genomic Footprinting**: Interactive Tn5 transposase cleavage profile curves showing TF binding protection.
  - **GRN Benchmark & Evaluation**: Precision-Recall curves & confusion matrix evaluating inferred GENIE3 edges vs ATAC peak intersections (with On-Demand evaluation trigger).
* **Selection & Controls**: TF / Motif search, cluster comparison dropdown, promoter window ($\pm 2\,\text{kb}, \pm 10\,\text{kb}$).

---

#### 🪟 Supporting Context: Docked Inspector & Quick-Checker

To avoid cluttering the regulatory tabs, non-regulatory inspection lives in a persistent, collapsible **Context Drawer**:
- **Cell Inspector**: Live QC metrics (counts, genes, mt%), cell type / Leiden ID, live view value.
- **Cluster Inspector**: Cell type composition (Tangram deconvolution), cluster marker preview, centrality scores.
- **Gene Expression & DGEA Quick-Checker**: Instant single-gene search & differential expression table/volcano plot without losing the active regulatory tab.

---

#### 🔗 Cross-Tab Interconnections

1. **Tab 1 $\longleftrightarrow$ Tab 4 (The Validation Loop)**:
   - Inspecting TF *FOXM1* in Tab 1 $\to$ Click *"Validate in ATAC"* $\to$ Opens Tab 4 with *FOXM1* motif accessibility and footprint pre-selected.
   - Dual-map view: **Chromatin Accessibility** (chromVAR) vs. **Transcriptional Activity** (AUCell) side-by-side.
2. **Tab 1 $\longleftrightarrow$ Context Drawer (Target Validation)**:
   - Clicking any target gene in a regulon instantly previews its spatial expression in the quick-checker.
3. **Tab 3 $\longrightarrow$ Tab 1 (Signaling to TF Activation)**:
   - Selecting a Ligand-Receptor axis in Tab 3 $\to$ Click *"Show Activated TFs"* filters Tab 1 to downstream transcription factors.
4. **Tab 2 $\longleftrightarrow$ Tab 1 (ceRNA vs. Transcription)**:
   - Clicking a ceRNA hub in Tab 2 cross-checks if it is regulated by a TF in Tab 1.

---

#### 🛡️ Modality-Aware Handling (No Ghost Tabs)

| Modality | When Present | When Absent (e.g. Visium / Xenium without ATAC / SPONGE) |
| :--- | :--- | :--- |
| **GENIE3 (Tab 1)** | Fully active with TF leaderboard, targets, and AUCell maps. | Informative card: *"No GENIE3 GRN in dataset. Upload network to enable."* |
| **SPONGE (Tab 2)** | Fully active with ceRNA modules & $mscor$ matrix. | Informative card: *"ceRNA network not present. Upload SPONGE network to enable."* |
| **LIANA+ (Tab 3)** | Spatial LR pairs, cosine similarities, and co-occurrence. | Spatial neighbor statistics (Squidpy co-occurrence) remain active. |
| **ATAC (Tab 4)** | chromVAR, Differential Motifs, Footprints, GRN Evaluation. | Clean banner: *"Epigenetic validation requires Multiome (paired RNA+ATAC). Connect Multiome data to unlock."* |


## Color by UX Decision 

Currently, the hexmap coloring can be invoked either by with the dropdown or with a tab. This is redundant. The only thing currently only in the dropdown are the "other" obs columns, which could easily be added in the Cell Information tab. 

### Option A: Keep both
Not a fan though 

### Option B: Keep only tabs
Simply remove the dropdown. 

### Option C: Un-synchronize the coloring and the tabs
Maybe there is a use case for coloring the map by X but looking at the table for Y. However this might be better solved by the compare view. 

## Evaluation of UX Options for Cell vs. Cluster Selection

Context:  
When interacting with the spatial transcriptomics map (hexagons / spots), users click on individual hexagon units to inspect biological properties, cluster assignments, quality control metrics, and regulatory network activities.

Currently, clicking behavior is modal (dependent on `selectedView`):
- When the map is colored by **Leiden** (`selectedView === 'leiden'`), clicking a hexagon selects and outlines the **entire Leiden cluster** and jumps to the **Cluster Information** tab.
- When the map is colored by **any other property** (e.g. `cell_type`, a gene, or a regulatory score), clicking a hexagon selects only that **individual cell** and jumps to the **Cell Information** tab.

---


### Option A: Cell-First with Direct 1-Click Cluster Navigation *(Recommended)*
- **Behavior**: Clicking any hexagon **always** selects that specific cell and opens **Cell Information**.
- **Cluster Integration**: Inside the Cell Information overview card, the `Cluster [X]` badge is rendered as an interactive button/chip (`[🏷️ Cluster 3 →]`). Clicking it highlights the cluster and switches to **Cluster Information**.
- **Pros**:
  - Predictable and consistent: a cell click always shows that cell's information regardless of which color layer is displayed.
  - Users can easily inspect an individual cell even when the map is colored by Leiden.
  - Preserves 1-click access to the full cluster view.
- **Cons**:
  - Requires one extra click to outline the whole cluster if the user only wanted cluster-level inspection.

---

### Option B: Unified Selection (Simultaneous Cell + Cluster Selection)
- **Behavior**: Clicking a cell selects **both** the specific cell and its parent cluster.
  - The clicked cell receives a solid border highlight.
  - The surrounding cluster members receive a subtle or dashed outline.
- **Tab Behavior**: Tab switching is non-aggressive (the user stays on whatever tab they are currently viewing, e.g. Differential Motif Activity, Cell Info, or Cluster Info, and all tabs update reactively with the new selection).
- **Pros**:
  - Context is never lost: both single-cell properties and cluster-level stats are available simultaneously.
  - Prevents disorienting automatic tab jumps.
- **Cons**:
  - Visual density on the map is higher (showing both single-cell outline and cluster outline).

---

### Option C: Explicit Mode Toggle (`Cell Mode` vs. `Cluster Mode`)
- **Behavior**: A toggle button group on the map toolbar: `🔘 Cell Selection` | `🔘 Cluster Selection`.
  - In **Cell Mode**: clicks highlight individual cells and open Cell Information.
  - In **Cluster Mode**: clicks outline all cells in the cluster and open Cluster Information.
- **Pros**:
  - Explicit user intent with no ambiguity.
- **Cons**:
  - Extra UI control; user must remember which mode is active.

---

## Dataset Generalization Rules

1. **Leiden Clustering (`leiden`)**:
   - Always computed during preprocessing / upload for all datasets.
   - Guaranteed primary categorical fallback when `cell_type` is absent.
   - Cluster IDs are normalized to integer numbers for strict equality checking.

2. **Cell Type Annotation (`cell_type`)**:
   - Optional (present in annotated datasets like Visium reference sets, absent in raw or unannotated single-cell / Xenium uploads).
   - Fallback hierarchy for default map coloring:
     1. `cell_type` (if present in `obs`)
     2. `leiden` (always present in `obs`)
     3. `regulatory_scores` (if computed)
     4. First numeric or categorical feature in `obs`
   - UI elements (e.g. Cell Type Distribution table in Cluster Information) gracefully handle datasets without `cell_type` by displaying "Unannotated" or hiding the distribution card.

3. **Cell Overview Card**:
   - Title defaults to `cell_type` if available, otherwise `Cluster <leiden_id>` or `Cell <barcode_short>`.
   - Cluster chip is always rendered whenever `leiden` is defined.

---

## Which AnnData should scores be computed on? (Tangram `if`/`else`)

**Status: open — needs a decision.** Documented behaviour and evidence in
[DEVELOPER.md](DEVELOPER.md) under *"Tangram runs score the projected object, not the measured
one"*. This entry is only the decision itself.

### The current behaviour

[calc_scores.py:371-456](backend/calc_python_scores/calc_scores.py#L371-L456) scores exactly one
object, chosen by an `if`/`else`:

```python
if args.tangram:
    ... compute_spatial_scores(adata_tangram, "tg", ...)   # Tangram-projected
else:
    compute_spatial_scores(adata_work, "st", ...)          # measured spatial
```

Multiome forces `use_tangram = True`, so **every multiome dataset is scored only on projected
expression**. This is intended, not a slip — FoPra p.38 states it explicitly: *"If Tangram was
enabled, these scores were calculated on the Tangram-derived AnnData object. Otherwise, they
were calculated directly on the [spatial] object."*

### Why it needs revisiting

1. **The data model asks for a file the design never produces.**
   `MultiomeDataset.validate_paths()` lists `adata_st_scores` as `always_required`, but multiome
   implies Tangram, so `"st"` is never computed. On the heart builtin this was "resolved" by
   hand-copying the projected object to `adata_st_scores.h5ad` — byte-identical, verified by
   md5. Anything reading that filename believed it had measured data.

2. **Two resolution functions disagree about which file to use.**

   | Function | Used for | Multiome picks |
   |---|---|---|
   | `_determine_adata_path` (dataset_management.py) | registration, GeoJSON source | `adata_tg_scores_path` |
   | `_resolve_adata_path` (main.py) | live API fetches (`/api/obsm_tables`, gene values) | `adata_st_scores_path` |

   So the map's embedded properties and the live-fetched values can come from *different files*.
   The byte-identical copy hid this; making the two files honest surfaced it immediately as
   "Regulatory Scores greyed out". Whatever is decided below, **these two must agree.**

3. **The projected object is not preprocessed for clustering or autocorrelation.**
   `project_genes` is never followed by re-normalisation: `expm1(X)` row sums have CV 12.9
   (measured spatial: exactly 10000, CV 0.0000), there is no `log1p` record and no `counts`
   layer, and values reach 28.97 so `expm1` ≈ 3.8e12 — not log-scale. Library size therefore
   re-enters as a leading PCA component. Separately, projection is spatial smoothing, which
   systematically **inflates** Moran's I / Geary's C, so autocorrelation on projected values is
   not comparable to autocorrelation on measurements.

4. **It is not user-visible.** Nothing in the UI says the Gene Expression tab is showing imputed
   rather than measured values.

### Options

**A. Score both objects (drop the `else`).** One-line change; `validate_paths()` becomes
honest; enables a real measured-vs-projected comparison. Costs roughly double the score runtime
for Tangram runs and doubles output size. Still needs a decision on which one the UI defaults
to, and the two resolution functions must be aligned.

**B. Keep scoring only one, but make it explicit and consistent.** Align
`_resolve_adata_path` with `_determine_adata_path`, drop `adata_st_scores` from
`always_required` for multiome, and label the source in the UI (e.g. "values are
Tangram-imputed"). Cheapest option; keeps the current science, fixes the incoherence and the
mislabelling.

**C. Score the measured object and project only what needs projecting.** chromVAR genuinely
requires the mapping (it is single-cell ATAC → spots). Expression-derived scores do not — they
could run on measurements, with Tangram used solely for cell-type composition and chromatin
projection. Most defensible scientifically; largest change.

**D. Let the user choose per upload.** A "compute scores on: measured / Tangram-projected /
both" option. Most flexible, most surface area, and multiplies the support burden.

### Notes for whoever picks

- The "Tangram mapped datasets" dropdown category exists but is currently a **no-op**: it swaps
  `adata_path` for `tangram_adata_path`, which are the same file for multiome. It also cannot
  work as-is, because `geojson_path` is spread through unchanged (one GeoJSON per dataset id)
  and the backend resolves paths from the registered dataset by `dataset_id`, ignoring the
  frontend's swapped field. A genuine toggle needs the two variants registered as separate
  dataset ids.
- Option A or C would make that dropdown meaningful for the first time.

### Correction to "Dataset Generalization Rules" above

That section states Leiden is *"always computed during preprocessing / upload for all
datasets"*. It is not. `clustering()` runs only when a Squidpy score that needs a cluster key
was requested **and** `leiden` is absent **and** the relevant cluster-key argument is literally
`"leiden"` ([calc_scores.py:91-96](backend/calc_python_scores/calc_scores.py#L91-L96)). The
heart builtin shipped with `scores.squidpy = false` and therefore had no `leiden` at all, which
is why its Cluster Information and Co-occurrence tabs were empty. Any UI logic that assumes
`leiden` is always present needs a fallback.
