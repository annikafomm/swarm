# SWARM Design Decisions & UX Options

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
