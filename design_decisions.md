# SWARM Design Decisions & UX Options

## Core Philosophy: Regulatory-First Architecture

**Problem**: The tool currently has 14 flat, unstructured tabs. Half of them are disabled or empty on datasets lacking ATAC or SPONGE, and generic features (QC, cell info, single-gene expression) overshadow SWARM's core novelty.

**Guiding Principle**: SWARM's primary scientific novelty is **mapping, scoring, and validating inferred gene regulatory networks (GENIE3), ceRNA crosstalk (SPONGE), and chromatin states (ATAC) on spatial coordinates**. Standard spatial tasks (cell typing, QC, single-gene inspection) serve as **supporting context**, not the hero experience.

---

## The 4 Regulatory Pillar Tabs + Context Drawer

Instead of 14 flat tool tabs, the interface is organized into **4 core regulatory question tabs** with a persistent **Context & Inspector** drawer:

```
┌────────────────────────────────────────────────────────────────────────────────────────────────────────┐
│  [ 🌲 1. Transcriptional GRNs ]  [ 🧽 2. ceRNA Networks ]  [ 🔄 3. Niche & Signaling ]  [ 🔓 4. Epigenetic Validation ]  │
└────────────────────────────────────────────────────────────────────────────────────────────────────────┘
```

### 🌲 Tab 1: Transcriptional GRNs (GENIE3 & Regulons)
* **Goal**: Explore which transcription factors drive spatial gene programs and discover their downstream targets.
* **Contents**:
  - **Regulon Leaderboard**: Ranked TFs $\to$ 1-click project AUCell / ULM regulon score onto the spatial map.
  - **Target Gene Explorer**: Edge-weighted target gene list ($w_{ij}$) for the active TF.
  - **Cell Composition $\times$ TF Activity**: Sender cell type abundance $\to$ receiver TF activation.
* **Selection & Controls**: Score metric (`AUCell`, `Moran's I`, `Geary's C`, `ULM`), TF search, min regulon size filter, color scale palette & percentile clipping.

### 🧽 Tab 2: Post-Transcriptional ceRNA Networks (SPONGE)
* **Goal**: Explore microRNA-mediated competitive endogenous RNA crosstalk and post-transcriptional spatial modules.
* **Contents**:
  - **Spatial ceRNA Module Activity**: spongeffects scores (`AUCell`, `GSVA`, `ssGSEA`) on spatial spots.
  - **ceRNA Interaction Matrix**: Significant gene-gene pairs ($mscor$, partial correlation, $p$-adj).
  - **Hub ceRNAs**: Most central sponging RNAs across the tissue.
* **Selection & Controls**: Scoring algorithm, $mscor$ cutoff slider, FDR $p$-value threshold.

### 🔄 Tab 3: Spatial Signaling & Niches (LIANA+)
* **Goal**: Identify intercellular communication that triggers intracellular regulatory programs.
* **Contents**:
  - **Spatial Ligand-Receptor Table**: Ranked by Cosine Similarity, $p$-value, categories (e.g. *High L - High R*), global LR Moran's I.
  - **Spatial Co-occurrence & Neighborhoods (Squidpy)**: Co-occurrence probability curves over distance radius $r$, neighborhood enrichment heatmap.
* **Selection & Controls**: Source $\to$ Target cell type filters, category filter, distance scale slider ($r$).

### 🔓 Tab 4: Epigenetic & Benchmark Validation (Multiome ATAC)
* **Goal**: Validate whether inferred GRNs reflect real, physical DNA-binding and accessible chromatin.
* **Contents**:
  - **chromVAR Motif Accessibility**: Spatial deviation z-scores (JASPAR 2024 CORE) + Moran's I / Geary's C.
  - **Differential Motif Activity**: Cluster-specific enriched motifs ($\Delta\text{accessibility}$, FDR).
  - **Genomic Footprinting**: Interactive Tn5 transposase cleavage profile curves showing TF binding protection.
  - **GRN Benchmark & Evaluation**: Precision-Recall curves & confusion matrix evaluating inferred GENIE3 edges vs ATAC peak intersections (with On-Demand evaluation trigger).
* **Selection & Controls**: TF / Motif search, cluster comparison dropdown, promoter window ($\pm 2\,\text{kb}, \pm 10\,\text{kb}$).

---

## 🪟 Supporting Context: Docked Inspector & Quick-Checker

To avoid cluttering the regulatory tabs, non-regulatory inspection lives in a persistent, collapsible **Context Drawer**:
- **Cell Inspector**: Live QC metrics (counts, genes, mt%), cell type / Leiden ID, live view value.
- **Cluster Inspector**: Cell type composition (Tangram deconvolution), cluster marker preview, centrality scores.
- **Gene Expression & DGEA Quick-Checker**: Instant single-gene search & differential expression table/volcano plot without losing the active regulatory tab.

---

## 🔗 Cross-Tab Interconnections

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

## 🛡️ Modality-Aware Handling (No Ghost Tabs)

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
