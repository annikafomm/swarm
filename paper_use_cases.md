# Paper use cases — heart (AVN multiome) and mouse visual cortex

Research output. Grounded in the two builtin datasets as they actually exist on disk (verified,
not assumed), a literature sweep per tissue, and a methods-landscape review. Dataset facts and
pipeline caveats are documented in [DEVELOPER.md](DEVELOPER.md); the open pipeline question is
in [design_decisions.md](design_decisions.md).

---

## Recommendation up front

**Write a methods/tool paper, not a biology paper, and use the two datasets as complementary
demonstrations rather than as two results chapters.**

Both datasets are **n=1** — one donor/one section each. That is fatal for any biological claim
that needs to generalise, and no amount of spots fixes it: adjacent Visium spots are spatially
autocorrelated by construction, so treating 1859 vs 115 spots as sample sizes is
pseudoreplication, not power. Every biological statement has to be phrased "in this donor".

What the datasets *can* carry is a methods paper whose contribution is **honesty about
regulatory scoring in space**: the same regulon set, scored four ways, on two tissues, with the
disagreements and the provenance of each score surfaced as first-class results. Three of the
strongest candidate contributions below are diagnostics, not discoveries — and that is the
paper's best available shape.

The two datasets are unusually well matched for that: human vs mouse, multiome vs
transcriptome-only, expert anatomical annotation vs deconvolved layers, and one tissue where
unsupervised domain detection recovers anatomy and one where it demonstrably does not.

---

## Killed before proposing — do not claim these

The literature sweep refuted four things, including one I had assumed was the standout angle.

1. **"Projecting chromVAR motif activity into space is novel."** Refuted by the paper whose
   method we use: Biancalani et al. 2021 (*Nature Methods*, Tangram) already maps SHARE-seq to
   space in Fig. 5.
2. **"First tool combining spatial regulatory-network inference with multiome validation."**
   Dead. ISON (*Nature Communications*, June 2026) does exactly this; STARNet (2025), SpaGRN
   (*Cell Systems* 2025), SCRIPro and SCENIC+ occupy adjacent ground.
3. **Circular validation in the heart.** chromVAR motif activity *and* the RNA-derived
   TF/regulon scores are placed in space by the **same** Tangram mapping. Their spatial
   co-localisation is therefore partly an artefact of shared projection, not independent
   corroboration. Any "chromatin confirms expression" claim must bound this explicitly.
4. **IEG dominance in visual cortex is not a discovery.** I had flagged EGR1/FOS/JUNB/NR4A1 as
   the headline. But the GENIE3 network for this dataset was learned from **GSE216973, a bulk
   light/dark ZT time-course of mouse V1** (Apelblat et al.). The regulon set is IEG-heavy
   *because the training design was*. Presenting it as a biological finding would be circular.
   It remains usable as a *substrate* for methodological questions, not as a result.

---

## Ranked use cases

### PRIMARY 1 — Estimator, not tissue, decides which regulons look "spatially structured"

*Both datasets. Methodological.*

**Question.** Is spatial autocorrelation of regulon activity a transferable readout, or does the
arbitrary choice of enrichment estimator decide the answer?

**Why this is the lead.** It is the one contribution already *empirically confirmed* on the real
data. Scoring identical GENIE3 regulons on identical coordinates, per-spot variance differs by
roughly two orders of magnitude between estimators:

| estimator | mean SD across regulons (heart) | value range |
|---|---|---|
| GSVA | **0.278** | −0.82 → 0.92 |
| VIPER | 0.048 | −4.75 → 10.5 |
| ssGSEA | 0.005 | −0.25 → 0.75 |
| AUCell (GENIE3) | 0.0028 | 0.34 → 0.85 |
| AUCell (SPONGE) | **0.0013** | 0.73 → 0.92 |

No column is constant (checked: 0 zero-variance columns), so this is not a defect — AUCell and
ssGSEA are simply near-flat spatially on this data where GSVA is not. Any tool that silently
picks one estimator is silently picking its own conclusions. **Every** spatial regulon tool
(SpaGRN, SCRIPro, STARNet, ISON, SCENIC+) picks one silently; decoupler benchmarked estimators
against perturbation ground truth but never against *spatial* rankings.

**Plan.** Rank the 20 regulons by Moran's I under each of the four estimators, in both tissues;
all 6 pairwise Spearman correlations of those rankings; regress each estimator's I on regulon
size, mean expression and detection rate to identify the mechanism of disagreement. Replace the
shipped analytic null — the pipeline runs `spatial_autocorr` with `n_perms=None`, so only a
normal approximation exists — with a size- and expression-matched permutation null.

**Figure.** (a) four-estimator rank concordance heatmap per tissue; (b) scatter of I under
AUCell vs GSVA with the regulons labelled; (c) variance-per-estimator bar chart (the table
above); (d) mechanism regression; (e) same two panels for the measured-vs-projected heart layers.

**Caveats.** GSVA and ssGSEA both come through the same spongEffects wrapper, so the "four
estimators" are really ~3 independent families. VIPER covers 19 of 20 regulons, not all.

---

### PRIMARY 2 — How much of the spatial regulatory signal is Tangram?

*Heart. Methodological, with a within-section control no other tool ships.*

**Question.** How much of the apparent spatial structure of regulatory scores is tissue
organisation versus mapping-induced smoothing?

**Why it matters.** The imputation literature establishes that imputation inflates Moran's I,
but nobody has quantified it *for regulatory scores* — regulons, ceRNA modules, motif activity —
and no spatial GRN tool ships this control. This dataset is unusually able to answer it because
**both layers exist for the same 3153 spots**: the Tangram-projected object (26,868 genes) and
the measured spatial object (33,538 genes, log1p, raw counts preserved).

**Plan.** Align on shared barcodes and genes; per-gene and per-spot Spearman between layers;
ΔMoran's I per gene; re-run node-vs-rest and cardiac_skeleton-vs-fibrosis DGEA on both layers
and report rank overlap; recompute all 20 regulons and 16 ceRNA modules identically on both.
Then isolate the smoothing component directly: shuffle spot labels within `adata_map.h5ad`,
re-project, recompute module-level Moran's I, and report the inflation attributable to the
mapping alone.

**Caveats.** The two layers use different normalisation, so scores must be *recomputed*
identically on both, never compared as-shipped. Gene-identifier namespaces differ between the
layers (Ensembl vs symbol) and must be reconciled before any join, or every concordance number
is wrong. SPONGE modules are named by Ensembl ID and are biologically uninterpretable until
resolved.

---

### PRIMARY 3 — Where does the node end? Boundary-resolved regulatory gradient across the AV insulating plane

*Heart. The strongest biological framing available.*

**Question.** Is the nodal regulatory programme bounded by a sharp step at the
node/cardiac_skeleton interface, or does it decay as a penumbra into peri-nodal myocardium — and
over how many spot widths?

**Why it is interesting.** This is the molecular form of the century-old "transitional cell"
question in the AV conduction axis. Kanemaru 2023 reported cell-state composition per
expert-annotated AVN structure and built pySCENIC GRNs, but never parameterised regulatory
activity *by distance to the boundary*; all TBX3/TBX5 nodal-identity work is mouse,
genotype-driven and non-spatial. Recasting niche-vs-rest as a distance-parameterised gradient is
the novel move, and it is what an interface-centric tool is for.

**Plan.** Build a boundary covariate from `obsp['spatial_connectivities']` — graph steps to the
nearest node/cardiac_skeleton interface. Pre-specify the readout panel from the conduction-system
literature *and check motif availability against the actual 719-motif panel before testing*. For
each readout, compare a step model against a monotone-decay model by AIC. Generate the null by
**rotation/torus permutation of the label field and block permutation of contiguous patches** —
not spot-wise Wilcoxon.

**Caveats — and one is close to fatal.** Visium spot width (55 µm, 100 µm pitch) is comparable
to the width of the compact node itself, so step-vs-gradient may be **formally unresolvable**;
the paper needs a panel that says so rather than forcing an answer. And Tangram projection is
spatial smoothing, so it *manufactures* exactly the graded appearance the analysis looks for —
this use case is only publishable *alongside* PRIMARY 2 as its control.

---

### SUPPORTING 4 — Two collagen compartments: is the annulus fibrosus regulatorily separable from reactive fibrosis?

*Heart. Runs a contrast the pipeline never ran.*

Do the developmentally programmed cardiac skeleton (719 spots) and reactive fibrosis (138 spots)
run distinct regulatory programmes, or the same collagen programme at different intensity? This
is the mechanism underneath age-related AV block. Kanemaru characterised FB4_activated
fibroblasts and explicitly noted they were *not* abundant in the AVN node structure, leaving
this open by design.

Primary test is effect-vector geometry: compute cardiac_skeleton-vs-myocardium and
fibrosis-vs-myocardium effect vectors over the 738 ULM TFs and compare their direction, rather
than testing each niche against "rest". Mandatory composition control using the 12
`tangram_ct_pred` fractions, since a fibroblast-abundance difference would mimic a regulatory
difference. Anchor discriminating TFs in chromatin via the already-computed
`collagen_fibril_orga` GRN evaluation.

**Caveats.** The "fibrosis" label may be annulus captured off-plane; there is no H&E
adjudication available. The candidate panel is asymmetric across layers (no SCX or RUNX1 motif,
no OSR1 in ULM), so "developmental vs reactive" cannot be tested symmetrically.

---

### SUPPORTING 5 — Continuous cortical depth, not discrete layers, for regulon activity

*Visual cortex. Best available biological framing there.*

L2/3 glutamatergic neurons in V1 form an experience-dependent transcriptomic *continuum*
covarying with pial depth (Cheng/Butrus/Zipursky). That result is about gene expression and
cell-type identity, established by MERFISH on ~500 targeted genes — with no regulon-level or
spatial-autocorrelation analysis. SCENIC+ reported *discrete* layer-specific Visium regulons.
Nobody has asked whether regulon **activity** is continuous or stepwise in depth.

Derive a per-spot pial-to-white-matter coordinate from the CLR-transformed laminar proportions,
fit a GAM of regulon score on depth, and compare a smooth-depth model against a discrete-layer
factor per regulon. Contrast identity regulons (NFIB, RFX3, SMAD5) against activity regulons
(EGR1, FOS, JUNB, NR4A1).

**Caveats.** Section geometry is the binding constraint — measured laminar centroids show the cut
is not a clean perpendicular column, so a whole-section depth model would smear layers. The depth
coordinate derives from the same deconvolution used as covariates, risking soft circularity.
L5 IT / L5 NP / L5 PT is a projection-class distinction, not a depth one, and must not be ordered
along the axis.

**External validation available:** Gray et al. 2017 (*eLife*) provides layer-resolved ATAC motif
enrichment for adult mouse V1 — an independent modality to check layer assignments against,
which partly compensates for this dataset having no ATAC of its own.

---

### SUPPORTING 6 — The cluster key is a provenance dimension

*Both. Turns a negative result into a stated precondition.*

Every cluster-keyed readout (centrality, neighbourhood enrichment, co-occurrence, differential
motif activity, DGEA) is conditional on a partition the user chooses. In the heart, unsupervised
Leiden **does not** recover the expert niches — ARI declines monotonically with resolution
(0.395 at k=2 down to 0.093 at k=12), and even at exactly 6 clusters ARI is 0.194, worse than
the 2-cluster solution. On its own that reads as a tool-choice error. Paired with a tissue where
the same pipeline and graph **do** recover anatomy, it becomes a precondition for the whole
spatial-domain literature.

Requires adding spatially-aware baselines (BayesSpace, STAGATE, GraphST) on the *same* Delaunay
graph, or a reviewer will say the negative result is just Leiden being the wrong tool.

---

### STRETCH 7 — Four-layer TBX3:TBX5 map of the human AV junction

No published human dataset displays TBX3/TBX5 antagonism spatially across the AV junction, and
the dataset holds four regulatory layers per spot (TF expression, chromVAR motif accessibility,
ULM activity, GENIE3 regulon score). Attractive — but **T-box motif family collinearity means
chromVAR cannot separate TBX3 from TBX5**, TBX3 has no GENIE3 regulon so the four-layer story is
complete only for TBX5, and regulon Moran's I values in this section are tiny and mostly
non-significant. Stretch because the discordance panel may be the only reportable outcome.

---

### STRETCH 8 — Geary's C adds nothing to Moran's I for regulon activity

The assumed "scale diagnostic" (high I + C≈1 = smooth gradient; low I + C≪1 = fine patchwork)
may be **degenerate by construction** for regulon scores, while remaining informative at gene
level on the identical graph. A clean, reusable rule if it holds. Risk: a reviewer dismisses it
as a known algebraic property of the two statistics for near-Gaussian row-standardised data.

---

## Rejected

- **Anything using the 16-spot vessel niche.** Drop it entirely rather than reporting with
  caveats — including it invites a reviewer to discount the whole differential-motif section.
- **Footprint-based claims.** Footprinting was requested in the original run but produced no
  outputs.
- **Cross-tissue quantitative comparison of score magnitudes.** Different species, network
  provenance, normalisation, and the projected-vs-measured issue make it invalid. Cross-dataset
  claims must be restricted to *within-dataset rankings*.
- **"Chromatin independently confirms expression" in the heart.** Circular — same Tangram
  mapping places both.
- **IEG biology as a finding in V1.** Circular — the network was trained on an IEG-inducing
  time-course.

---

## Statistical requirements (non-negotiable)

1. **Spots are not replicates.** Use rotation/torus permutation of label fields and block
   permutation of contiguous patches. Spot-wise Wilcoxon/t-tests are pseudoreplication.
2. **Replace the analytic null.** The pipeline computes autocorrelation with `n_perms=None`, so
   only a normal approximation exists. Build size- and expression-matched permutation nulls.
3. **n=1 phrasing throughout.** "In this donor" / "in this section". Permutation p-values test
   spatial randomness *within* a section, not reproducibility across animals.
4. **Report effect sizes first**, p-values second, given the pseudoreplication exposure.

## Gaps worth closing, by ceiling raised

| Fix | Unlocks |
|---|---|
| A second heart donor / section | Any generalisable AVN biological claim; converts PRIMARY 3 from descriptive to inferential |
| Resolve the projected-vs-measured pipeline question | Removes the confound under every heart autocorrelation number |
| ATAC for visual cortex | Separates TF expression from TF activity; currently impossible there |
| Regenerate V1 from uncorrupted mm10 symbols | **Done** — `scripts/restore_var_symbols.py` copied symbols positionally from `visual_cortex_visium_1_prepro_final.h5ad` (same order, no upper-case collisions), repairing `var_names` plus the gene-keyed `uns` score tables. **Do not fix by re-capitalising** — `str.upper()` is not invertible: correct mouse casing follows no single rule, and neither `.title()` nor `.capitalize()` recovers `Rb1cc1`, `4732440D04Rik`, `mt-Nd1` and `AI597479` together. Root cause fixed upstream: LIANA now uses the `mouseconsensus` resource instead of upper-casing the data |
| A GENIE3 network for V1 from a non-IEG-inducing design | De-circularises any activity-dependent claim |
| Resolve SPONGE module names to symbols | Makes the ceRNA modules biologically interpretable |

## Reviewer objections to prepare for

| Objection | Response |
|---|---|
| "n=1, so nothing is shown." | Concede for biology; the contributions are diagnostics validated within-section, plus a cross-tissue replication of the *methodological* finding. |
| "Spatial chromVAR projection isn't new." | Agreed — cite Biancalani Fig. 5. The claim is the *audit* of what projection does to autocorrelation, not the projection. |
| "ISON/STARNet already do spatial GRN + multiome." | Agreed. Differentiator is the estimator-concordance and provenance diagnostics, and per-edge falsification vs their aggregate AUROC. |
| "Chromatin/expression agreement is circular." | Stated explicitly as a bound, with the shuffled-mapping null quantifying it. |
| "Your negative clustering result is just bad Leiden." | Spatially-aware baselines on the same graph, plus the contrasting tissue where the pipeline succeeds. |
| "Estimator disagreement is a known artefact." | Known in the abstract; unquantified for *spatial* rankings, and the ~100× variance spread here is the concrete demonstration. |
