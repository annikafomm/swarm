# like Noah's run_liana.py file, but slight structure changes

import decoupler as dc
import pandas as pd
import scanpy as sc
from liana_folder.liana_bivariate import (
    cell_comp_tf_activity_similarity,
    ligand_receptor_relationships,
)
from liana_folder.liana_unsupervised import pathway_activities

def looks_like_mouse(adata, sample_size=500):
    """Guess the organism from gene-symbol casing.

    Mouse and human nomenclature differ by convention: mouse is `Nfib`, human is `NFIB`.
    So a var index that is overwhelmingly not upper-case is mouse (or another
    non-human organism using the same convention).

    Deliberately conservative — it only reports mouse when the evidence is strong, because
    the cost of guessing wrong is running the wrong ligand-receptor resource. Callers can
    override with an explicit `organism=`.
    """
    names = [str(n) for n in adata.var_names[:sample_size] if str(n)]
    if not names:
        return False
    lower_case_present = sum(1 for n in names if n != n.upper())
    return lower_case_present / len(names) > 0.5


def run_liana(
    adata,
    grn=None,
    pathway_net=None,
    cell_comp_key="tangram_ct_pred",
    organism=None,
):
    """Run the three LIANA+ analyses.

    `organism` is "human", "mouse", or None to auto-detect from symbol casing.

    Note on gene symbols: this function used to do
        adata.var.index = adata.var.index.str.upper()
    in place, so that the human-symbol `consensus` ligand-receptor resource would match. That
    mutation was permanent — the caller writes the same object back to disk — and it silently
    destroyed mouse symbols (`Xkr4` -> `XKR4`, `mt-Nd1` -> `MT-ND1`), breaking every join
    against a mouse annotation and making the symbols look human. It is not reversible either,
    because correct mouse casing is not a single pattern (`.title()` and `.capitalize()` each
    fail on different classes: `Rb1cc1`, `4732440D04Rik`, `mt-Nd1`, `AI597479`).

    It was also unnecessary. LIANA ships a mouse-native resource, `mouseconsensus`, and
    decoupler's `collectri`/`progeny` both take `organism=`, so the right fix is to request
    each reference network for the organism at hand rather than mangle the data to match a
    human one. `liana.resource.get_hcop_orthologs()` + `translate_resource()` cover other
    organisms via HCOP if one is ever needed.

    Note that a caller-supplied `grn`/`pathway_net` CSV is used verbatim -- it is assumed to
    already match the data (the visual-cortex builtin's GENIE3 network, for instance, is
    correctly mouse-cased).
    """
    adata.uns["liana_columns"] = {}
    adata.var.index.name = None

    if organism is None:
        organism = "mouse" if looks_like_mouse(adata) else "human"
        print(f"[liana] organism auto-detected from symbol casing: {organism}")
    else:
        print(f"[liana] organism specified by caller: {organism}")

    # Every reference network here is keyed by gene symbol, so all three must be requested for
    # the right organism. Getting this wrong does not raise -- it just silently scores against
    # almost no overlapping genes. Measured on the mouse visual-cortex builtin (31,053 symbols):
    #   progeny(human)   12/17,610 target genes present (0.1%)   vs  mouse 15,738/16,706 (94.2%)
    #   collectri(human) 13/6,675   target genes present (0.2%)   vs  mouse  6,421/6,582  (97.6%)
    #   consensus         9/2,017   LR genes present     (0.4%)   vs  mouseconsensus 1,559/1,741 (89.5%)
    if grn is None:
        grn = dc.op.collectri(organism=organism)
    else:
        grn = pd.read_csv(grn)

    if pathway_net is None:
        pathway_net = dc.op.progeny(organism=organism)
    else:
        pathway_net = pd.read_csv(pathway_net)

    resource_name = "mouseconsensus" if organism == "mouse" else "consensus"
    print(f"[liana] ligand-receptor resource: {resource_name}")

    ligand_receptor_relationships(adata, resource_name=resource_name)
    if cell_comp_key in adata.obsm.keys():
        cell_comp_tf_activity_similarity(adata, grn, cell_comp_obsm_key=cell_comp_key)
    pathway_activities(adata, pathway_net)

    return adata
