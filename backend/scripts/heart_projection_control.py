#!/usr/bin/env python3
"""
Does the chromatin-expression spatial agreement in the heart survive breaking the shared
Tangram projection?

Why this exists
---------------
`paper_use_cases.md` killed the heart's headline claim on circularity grounds:

    chromVAR motif activity *and* the RNA-derived TF/regulon scores are placed in space by the
    **same** Tangram mapping. Their spatial co-localisation is therefore partly an artefact of
    shared projection, not independent corroboration.

That objection is fatal to "chromatin confirms expression" as stated, and PRIMARY 3 (the AV-node
boundary gradient -- the strongest biological framing available for this dataset) is explicitly
only publishable alongside a control that bounds it. This script is that control, and it is
possible only because the heart builtin ships **both** layers for the same 3153 spots:

  * `adata_tg_scores.h5ad`        Tangram-projected expression (26,868 genes) + chromVAR + graph
  * `adata_spatial_measured.h5ad` the actual Visium measurement (33,538 genes)

The measured layer never passes through the projection. So regulon activity recomputed on it is
placed in space by the microscope, not by Tangram, while chromVAR is still placed by Tangram:

    r_projected = corr( regulon activity from PROJECTED RNA , chromVAR motif activity )
    r_measured  = corr( regulon activity from MEASURED  RNA , chromVAR motif activity )

`r_projected` contains signal plus shared-projection artefact. `r_measured` keeps the signal but
drops the shared-linear-combination path, because the two sides no longer come from the same
mapping matrix applied to the same snRNA cells. The drop bounds how much of the
published-style agreement the mapping manufactured.

Read the drop, not the absolute values:
  * drop ~ 0 and r_measured > 0  -> agreement is not a projection artefact; PRIMARY 3 survives
  * drop large, r_measured ~ 0   -> the agreement WAS the projection; the diagnostic is the result

Caveat that must survive into any write-up: the placements are not fully independent even here.
Tangram learns its mapping by matching snRNA against the *spatial* RNA, so chromVAR's placement
is informed by the same measurement the measured-layer regulons come from. This bounds the
artefact rather than eliminating it -- it removes the strong dependency, not every path.

Also reported, since it is the same measurement once aligned: per-regulon Moran's I on both
layers using the SAME graph. Mapping-induced smoothing is expected to inflate spatial
autocorrelation; `delta_morans_I` quantifies that for regulatory scores specifically, which the
imputation literature establishes in general but not for regulons.

Reading strategy -- why this does not use scanpy
-----------------------------------------------
`/workspaces/swarm` is an **NFS mount**. `sc.read_h5ad()` pulls the entire object -- every obsm
score matrix, PCA/UMAP embeddings, tissue images -- when this analysis needs perhaps 1% of it.
On a saturated mount that stalls indefinitely in uninterruptible `D` state (observed: a 128 MB
file not read in ten minutes, load average 232). So every read here is scoped: only the
regulon-target gene columns, the spatial graph, the chromVAR matrix, and the genesets dict.

`--stage` copies the inputs to local `/tmp` (overlay, not NFS) first. The copy still pays one
sequential NFS read, but every subsequent access is local -- worth it if you plan to re-run.

Usage
-----
    python3 heart_projection_control.py --inventory            # cheap: metadata only
    python3 heart_projection_control.py --out results.json
    python3 heart_projection_control.py --stage --out results.json
"""

import argparse
import json
import os
import shutil
import sys
import time

HEART = "/workspaces/swarm/backend/data/heart"
PROJECTED = f"{HEART}/plasmidpoop/adata_tg_scores.h5ad"
MEASURED = f"{HEART}/plasmidpoop/adata_spatial_measured.h5ad"
MOTIF_TO_TF = f"{HEART}/plasmidpoop/motif_to_tf.csv"
STAGE_DIR = "/tmp/swarm_projection_control"


def log(msg):
    print(f"[proj-control] {msg}", flush=True)


# ---------------------------------------------------------------------------
# scoped h5ad readers -- never load a whole object
# ---------------------------------------------------------------------------

def _dec(xs):
    return [x.decode() if isinstance(x, bytes) else str(x) for x in xs]


def _index_of(group):
    key = group.attrs.get("_index", "_index")
    if isinstance(key, bytes):
        key = key.decode()
    return _dec(group[key][:])


def _column(group, name):
    """One obs/var column, decoding AnnData's categorical encoding."""
    import h5py

    node = group[name]
    if isinstance(node, h5py.Group):
        cats = _dec(node["categories"][:])
        return [cats[c] if c >= 0 else None for c in node["codes"][:]]
    return _dec(node[:])


def _encoding(node):
    enc = node.attrs.get("encoding-type")
    if isinstance(enc, bytes):
        enc = enc.decode()
    if enc:
        return enc
    import h5py

    return "csr_matrix" if isinstance(node, h5py.Group) else "array"


def read_gene_columns(f, col_idx, n_obs):
    """Return a dense (n_obs, len(col_idx)) float64 array for the requested gene columns only.

    Handles the three encodings AnnData writes. The point is to transfer a few thousand columns
    instead of tens of thousands, which is what makes this runnable over NFS.
    """
    import numpy as np
    from scipy import sparse

    node = f["X"]
    enc = _encoding(node)
    order = np.argsort(col_idx)          # h5py needs increasing fancy indices
    sorted_idx = np.asarray(col_idx)[order]
    inverse = np.argsort(order)

    if enc == "array":
        # Dense: fancy-index the gene axis directly.
        block = node[:, sorted_idx]
        out = np.asarray(block, dtype=np.float64)

    elif enc == "csc_matrix":
        # Column-major: read only the requested columns' slices. Optimal case.
        indptr = node["indptr"][:]
        cols = []
        for c in sorted_idx:
            lo, hi = int(indptr[c]), int(indptr[c + 1])
            dense = np.zeros(n_obs, dtype=np.float64)
            if hi > lo:
                dense[node["indices"][lo:hi]] = node["data"][lo:hi]
            cols.append(dense)
        out = np.column_stack(cols) if cols else np.zeros((n_obs, 0))

    elif enc == "csr_matrix":
        # Row-major: columns are scattered, so read the sparse triple and slice in memory.
        # Still far cheaper than the whole file, since we skip obsm/obsp/uns/images entirely.
        shape = tuple(int(x) for x in node.attrs["shape"])
        M = sparse.csr_matrix(
            (node["data"][:], node["indices"][:], node["indptr"][:]), shape=shape
        )
        out = np.asarray(M[:, sorted_idx].todense(), dtype=np.float64)

    else:
        raise SystemExit(f"unsupported X encoding: {enc!r}")

    return out[:, inverse]


def read_sparse(node):
    """Read an obsp/obsm sparse matrix in full (these are small: spots x spots or spots x k)."""
    import numpy as np
    from scipy import sparse

    enc = _encoding(node)
    if enc == "array":
        return sparse.csr_matrix(np.asarray(node[:]))
    shape = tuple(int(x) for x in node.attrs["shape"])
    data = (node["data"][:], node["indices"][:], node["indptr"][:])
    M = sparse.csr_matrix(data, shape=shape) if enc == "csr_matrix" else \
        sparse.csc_matrix(data, shape=shape)
    return M.tocsr()


def read_obsm_frame(f, key):
    """Read one obsm entry as (array, column_names)."""
    import numpy as np
    import h5py

    if "obsm" not in f or key not in f["obsm"]:
        return None, None
    node = f["obsm"][key]
    if isinstance(node, h5py.Group) and node.attrs.get("_index") is not None:
        cols = node.attrs.get("column-order")
        cols = _dec(cols) if cols is not None else []
        arr = np.column_stack([np.asarray(node[c][:], dtype=np.float64) for c in cols])
        return arr, cols
    if isinstance(node, h5py.Group):
        return read_sparse(node).toarray().astype(np.float64), None
    return np.asarray(node[:], dtype=np.float64), None


def read_genesets(f, key="genie_genesets"):
    """Read uns['genie_genesets'] as {regulator: [targets]}."""
    import h5py

    if "uns" not in f or key not in f["uns"]:
        return {}
    node = f["uns"][key]
    if not isinstance(node, h5py.Group):
        return {}
    out = {}
    for reg in node:
        try:
            out[str(reg)] = _dec(node[reg][:])
        except Exception:
            continue
    return out


# ---------------------------------------------------------------------------
# analysis
# ---------------------------------------------------------------------------

def align_spots(proj_obs, meas_obs):
    """Match barcodes across layers, tolerating a sample-name prefix on one side."""
    def tail(b):
        return b.rsplit("_", 1)[-1]

    p = {tail(b): i for i, b in enumerate(proj_obs)}
    m = {tail(b): i for i, b in enumerate(meas_obs)}
    shared = sorted(set(p) & set(m))
    return [p[t] for t in shared], [m[t] for t in shared], shared


def zscore_mean(X):
    """Mean per-gene z-score across the gene axis -- the regulon score.

    Scale-free by construction, which is required here: the measured layer is library-size
    normalised to exactly 1e4 while the projected layer is a Tangram-weighted average of snRNA
    profiles that is NOT re-normalised (row sums vary ~13% CV). Any scale-sensitive score would
    compare normalisation rather than biology.
    """
    import numpy as np

    sd = X.std(axis=0)
    keep = sd > 0
    if keep.sum() < 2:
        return None, 0
    Z = (X[:, keep] - X[:, keep].mean(axis=0)) / sd[keep]
    return Z.mean(axis=1), int(keep.sum())


def rownorm_linear(X, label, target=1.0e4):
    """Rescale each spot linearly by its own row sum. No expm1 -- that is the whole point.

    Established empirically (see heart_projection_control_findings.md):

      * Tangram's projection is `X_space = M.T @ X_sc` (tangram/utils.py:368) -- a weighted SUM
        over reference cells, where M's rows sum to exactly 1.0. So the projected values are a
        LINEAR combination of the reference's already-log-normalised expression.
      * The projected layer's linear row-sum CV is 0.2194, matching the mapping matrix's
        per-spot column-sum CV of 0.2194 to four decimals. That per-spot factor is mapping mass:
        how many reference cells landed on the spot.
      * The measured layer's log-space row-sum CV is also 0.2194 -- Tangram maps more cells onto
        deeper spots, so it reproduces the real depth factor.

    Because the projected layer is a linear sum, that factor is MULTIPLICATIVE, so a per-gene
    z-score leaves it in every gene and every gene set inherits it -- which is what collapsed all
    20 regulons onto one field. Dividing it out linearly is the correct correction.

    Do NOT apply expm1 first. The projected values are sums of log-scale numbers, so expm1 is
    semantically meaningless, and numerically it is a disaster: it drives the row-sum CV from 0.22
    to 17.0 and leaves the single largest entry holding 99.15% of each row's mass (median),
    because expm1(29.6) is about 7e12.
    """
    import numpy as np

    totals = X.sum(axis=1, keepdims=True)
    cv = float(totals.std() / totals.mean()) if totals.mean() else float("nan")
    log(f"    {label}: linear row sums mean={totals.mean():.4g} CV={cv:.4f} -> rescaling")
    totals[totals == 0] = 1.0
    return X / totals * target


def normalize_like_counts(X, label):
    """Put a layer on a common per-spot scale: expm1 -> row-sum 1e4 -> log1p.

    This is the step my first attempt was missing, and its absence invalidated that attempt
    entirely. The measured layer is already log1p(CPM-1e4): its expm1 row sums are exactly 10000
    with CV 0.0000, so this transform is the identity for it. The Tangram-projected layer is NOT
    re-normalised -- its expm1 row sums average 2.1e10 with CV 12.86 -- so a per-gene z-score
    leaves an enormous per-spot factor in place, which every gene set then inherits. That is what
    collapsed all 20 regulons onto a single field (pairwise |r| 0.993, PC1 97.5%, effectively 1
    dimension out of 20).

    Normalising first removes that factor while keeping the estimator identical on both layers,
    which is what PRIMARY 2 requires: recompute scores identically, never compare as-shipped.
    """
    import numpy as np

    counts = np.expm1(np.clip(X, 0, None))
    totals = counts.sum(axis=1, keepdims=True)
    cv = float(totals.std() / totals.mean()) if totals.mean() else float("nan")
    log(f"    {label}: expm1 row sums mean={totals.mean():.3g} CV={cv:.4f}"
        + ("  (already normalised -- transform is ~identity)" if cv < 0.01
           else "  -> renormalising to 1e4"))
    totals[totals == 0] = 1.0
    return np.log1p(counts / totals * 1e4)


def aucell_scores(X, gene_names, genesets, tmin):
    """AUCell score per regulon, via decoupler, on a dense (spots x genes) matrix.

    Rank-based within each spot, so it is invariant to per-spot scaling -- which matters
    enormously here: the Tangram-projected layer is NOT library-normalised (expm1 row sums vary
    with CV 12.9, versus exactly 0.0 for the measured layer). A per-gene z-score removes gene
    scale but leaves that per-spot factor, so every gene set inherits the same dominant global
    field and all regulon scores collapse onto one another. AUCell does not have that failure
    mode, and it is also the estimator the SWARM pipeline itself uses.
    """
    import anndata as ad
    import decoupler as dc
    import numpy as np
    import pandas as pd

    net = pd.DataFrame(
        [(tf, g) for tf, ts in genesets.items() for g in set(ts) if g in set(gene_names)],
        columns=["source", "target"],
    )
    if net.empty:
        raise SystemExit("no regulon target present in the shared gene background")
    adata = ad.AnnData(
        X=np.asarray(X, dtype=np.float32),
        var=pd.DataFrame(index=pd.Index([str(g) for g in gene_names])),
    )
    dc.mt.aucell(adata, net=net, tmin=tmin, verbose=False)
    key = next((k for k in ("score_aucell", "aucell_estimate") if k in adata.obsm), None)
    if key is None:
        raise SystemExit(f"decoupler wrote no recognised AUCell key; obsm={list(adata.obsm)}")
    scores = adata.obsm[key]
    cols = list(scores.columns) if hasattr(scores, "columns") else []
    arr = scores.to_numpy() if hasattr(scores, "to_numpy") else np.asarray(scores)
    return {str(c): arr[:, i].astype(np.float64) for i, c in enumerate(cols)}


def morans_i(values, W):
    import numpy as np

    x = np.asarray(values, dtype=np.float64)
    z = x - x.mean()
    denom = float((z ** 2).sum())
    s0 = float(W.sum())
    if denom == 0 or s0 == 0:
        return float("nan")
    return (len(x) / s0) * (float(z @ (W @ z)) / denom)


def stage(paths):
    os.makedirs(STAGE_DIR, exist_ok=True)
    out = []
    for p in paths:
        dst = os.path.join(STAGE_DIR, os.path.basename(p))
        if os.path.exists(dst) and os.path.getsize(dst) == os.path.getsize(p):
            log(f"    already staged: {dst}")
        else:
            t = time.time()
            log(f"    staging {os.path.basename(p)} -> {dst} ...")
            shutil.copy2(p, dst)
            log(f"    staged in {time.time()-t:.0f}s")
        out.append(dst)
    return out


def inventory(projected, measured):
    import h5py

    for path, tag in ((projected, "PROJECTED (Tangram)"), (measured, "MEASURED (Visium)")):
        if not os.path.exists(path):
            log(f"MISSING: {path}")
            continue
        t = time.time()
        with h5py.File(path, "r") as f:
            obs_i, var_i = _index_of(f["obs"]), _index_of(f["var"])
            log(f"=== {tag} :: {os.path.basename(path)} ===")
            log(f"    {len(obs_i)} obs x {len(var_i)} var   (metadata read in {time.time()-t:.0f}s)")
            log(f"    obs e.g. {obs_i[:2]}")
            log(f"    var e.g. {var_i[:4]}")
            log(f"    var columns {sorted(f['var'])}")
            log(f"    X encoding  {_encoding(f['X'])}")
            log(f"    obsm        {sorted(f['obsm']) if 'obsm' in f else []}")
            log(f"    obsp        {sorted(f['obsp']) if 'obsp' in f else []}")
            log(f"    uns         {sorted(f['uns'])[:22] if 'uns' in f else []}")


def control(args, projected, measured):
    import h5py
    import numpy as np
    import pandas as pd
    from scipy import sparse, stats

    # ---- projected layer: metadata, genesets, graph, chromVAR ----
    t = time.time()
    log("opening projected layer (scoped reads only)")
    with h5py.File(projected, "r") as f:
        p_obs = _index_of(f["obs"])
        p_var = _index_of(f["var"])
        genesets = read_genesets(f)
        if "spatial_connectivities" not in f.get("obsp", {}):
            raise SystemExit("obsp['spatial_connectivities'] missing on the projected layer")
        W_all = read_sparse(f["obsp"]["spatial_connectivities"])
        chrom_arr, chrom_cols = read_obsm_frame(f, "chromvar_spot_scores")
        motif_names = _dec(f["uns"]["chromvar_motifs"][:]) if (
            "uns" in f and "chromvar_motifs" in f["uns"]) else []
        p_enc = _encoding(f["X"])
    log(f"    {len(p_obs)} obs x {len(p_var)} var, X={p_enc}, "
        f"{len(genesets)} regulons, graph {W_all.shape} ({W_all.nnz} edges) "
        f"in {time.time()-t:.0f}s")

    # projected-side alternative key: Ensembl IDs, if the object carries them
    with h5py.File(projected, "r") as f:
        p_ensembl = _column(f["var"], "ensembl_id") if "ensembl_id" in f["var"] else None

    # ---- measured layer: metadata only for now ----
    t = time.time()
    log("opening measured layer (scoped reads only)")
    with h5py.File(measured, "r") as f:
        m_obs = _index_of(f["obs"])
        m_var_index = _index_of(f["var"])
        m_var_cols = sorted(f["var"])
        m_candidates = {"index": m_var_index}
        for c in ("SYMBOL", "symbol", "gene_symbol", "feature_name", "gene_name"):
            if c in f["var"]:
                m_candidates[c] = _column(f["var"], c)
        m_enc = _encoding(f["X"])
    log(f"    {len(m_obs)} obs x {len(m_var_index)} var, X={m_enc}, "
        f"var columns {m_var_cols} in {time.time()-t:.0f}s")

    # ---- spots ----
    p_rows, m_rows, shared = align_spots(p_obs, m_obs)
    log(f"shared spots: {len(shared)} of {len(p_obs)} projected / {len(m_obs)} measured")
    if len(shared) < 100:
        raise SystemExit("too few shared spots -- barcode conventions do not match")

    # ---- gene namespace bridge ----
    # The two layers use different gene-identifier namespaces: the projected object is keyed by
    # HGNC symbol (with Ensembl in var['ensembl_id']), while the measured object is keyed by
    # Ensembl (with symbols in var['SYMBOL']). Either side can supply the join key, so try both
    # directions and take whichever overlaps more -- Ensembl is normally the safer of the two,
    # being stable and free of the symbol-casing hazards this dataset has already suffered.
    #
    # Regulons are named by SYMBOL regardless, so `needed` is resolved in symbol space and then
    # translated through whichever bridge wins.
    bridges = {}
    for name, vals in m_candidates.items():
        keys = [str(v) if v else "" for v in vals]
        bridges[f"symbol:measured.{name}"] = (list(p_var), keys)
    if p_ensembl is not None:
        bridges["ensembl:projected.ensembl_id"] = (
            [str(v) if v else "" for v in p_ensembl], [str(v) for v in m_var_index]
        )

    best, best_n, best_pair = None, -1, None
    for name, (pk, mk) in bridges.items():
        n = len({k for k in pk if k} & {k for k in mk if k})
        log(f"    bridge {name}: {n} shared keys")
        if n > best_n:
            best, best_n, best_pair = name, n, (pk, mk)
    log(f"    using bridge {best!r} ({best_n} shared keys)")
    p_keys, m_keys = best_pair

    # first occurrence wins, so a duplicated key maps to one column deterministically
    p_pos, m_pos = {}, {}
    for i, k in enumerate(p_keys):
        if k:
            p_pos.setdefault(k, i)
    for i, k in enumerate(m_keys):
        if k:
            m_pos.setdefault(k, i)
    shared_keys = set(p_pos) & set(m_pos)
    log(f"gene keys present in BOTH layers: {len(shared_keys)}")
    if len(shared_keys) < 1000:
        raise SystemExit("gene overlap too small -- namespace bridge failed")

    # Map symbol -> (projected col, measured col) so regulon definitions (symbol-named) resolve
    # regardless of which bridge won.
    sym_to_cols = {}
    for i, sym in enumerate(p_var):
        key = p_keys[i]
        if key and key in m_pos and sym not in sym_to_cols:
            sym_to_cols[sym] = (p_pos[key], m_pos[key])
    universe = set(sym_to_cols)
    log(f"regulon-addressable genes (symbol -> both layers): {len(universe)}")

    if not genesets:
        raise SystemExit("uns['genie_genesets'] missing -- cannot define regulons")

    # ---- which gene columns do we need? ----
    # For mean-z, only the regulon targets. For AUCell, the FULL shared gene set: AUCell scores a
    # gene set by where its members fall in each spot's ranking of ALL genes, so restricting the
    # background to regulon targets would distort every score.
    targets_only = sorted({g for ts in genesets.values() for g in ts if g in universe})
    if args.score == "aucell":
        needed = sorted(universe)
        log(f"AUCell needs the full shared gene background: {len(needed)} genes "
            f"({len(targets_only)} of them are regulon targets)")
    elif args.score in ("meanz-norm", "meanz-rownorm"):
        # Per-spot normalisation must use the whole shared background, or the row sums it
        # divides by would be computed over regulon targets only and would not be comparable
        # across layers with different gene coverage.
        needed = sorted(universe)
        log(f"{args.score} needs the full shared background for per-spot totals: "
            f"{len(needed)} genes ({len(targets_only)} are regulon targets)")
    else:
        needed = targets_only
        log(f"mean-z needs regulon targets only: {len(needed)} genes "
            f"({100*len(needed)/max(len(p_var),1):.1f}% of the projected gene axis)")
    if not needed:
        raise SystemExit("no regulon target gene is present in both layers")

    p_cols = [sym_to_cols[g][0] for g in needed]
    m_cols = [sym_to_cols[g][1] for g in needed]

    # ---- the only large reads, both scoped to `needed` ----
    t = time.time()
    log("reading projected expression for those genes only")
    with h5py.File(projected, "r") as f:
        Xp = read_gene_columns(f, p_cols, len(p_obs))
    log(f"    {Xp.shape} in {time.time()-t:.0f}s")

    t = time.time()
    log("reading measured expression for those genes only")
    with h5py.File(measured, "r") as f:
        Xm = read_gene_columns(f, m_cols, len(m_obs))
    log(f"    {Xm.shape} in {time.time()-t:.0f}s")

    # align rows to shared spots
    Xp = Xp[p_rows, :]
    Xm = Xm[m_rows, :]
    W = W_all[p_rows, :][:, p_rows]
    log(f"aligned: projected {Xp.shape}, measured {Xm.shape}, graph {W.shape}")

    col_of = {g: i for i, g in enumerate(needed)}

    # ---- chromVAR ----
    chrom_df = None
    if chrom_arr is not None:
        cols = chrom_cols or (motif_names if len(motif_names) == chrom_arr.shape[1]
                              else [f"motif_{i}" for i in range(chrom_arr.shape[1])])
        chrom_df = pd.DataFrame(chrom_arr[p_rows, :], columns=[str(c) for c in cols])
        log(f"chromVAR: {chrom_df.shape[1]} motifs over {chrom_df.shape[0]} shared spots")
    else:
        log("obsm['chromvar_spot_scores'] absent -- chromatin arm will be skipped")

    motif_tf = {}
    if os.path.exists(MOTIF_TO_TF):
        mt = pd.read_csv(MOTIF_TO_TF)
        log(f"motif_to_tf.csv columns: {list(mt.columns)}")
        mcol = mt.columns[0]
        tcol = next((c for c in mt.columns[1:] if mt[c].dtype == object), mt.columns[-1])
        for m, tf in zip(mt[mcol], mt[tcol]):
            if pd.isna(m) or pd.isna(tf):
                continue
            motif_tf.setdefault(str(tf).upper(), []).append(str(m))
        log(f"    {len(motif_tf)} TFs with at least one motif")
    else:
        log(f"motif_to_tf.csv absent at {MOTIF_TO_TF}")

    # ---- per-regulon ----
    # ---- score both layers with the SAME estimator ----
    if args.score == "meanz-rownorm":
        log("linear per-spot rescale of both layers (no expm1), then mean per-gene z")
        Xp = rownorm_linear(Xp, "projected")
        Xm = rownorm_linear(Xm, "measured ")
        proj_all = {tf: zscore_mean(Xp[:, [col_of[g] for g in ts if g in col_of]])[0]
                    for tf, ts in genesets.items()
                    if len([g for g in ts if g in col_of]) >= args.min_targets}
        meas_all = {tf: zscore_mean(Xm[:, [col_of[g] for g in ts if g in col_of]])[0]
                    for tf, ts in genesets.items()
                    if len([g for g in ts if g in col_of]) >= args.min_targets}
    elif args.score == "meanz-norm":
        log("normalising both layers to a common per-spot scale before scoring")
        Xp = normalize_like_counts(Xp, "projected")
        Xm = normalize_like_counts(Xm, "measured ")
        log("scoring with mean per-gene z on the normalised layers")
        proj_all = {tf: zscore_mean(Xp[:, [col_of[g] for g in ts if g in col_of]])[0]
                    for tf, ts in genesets.items()
                    if len([g for g in ts if g in col_of]) >= args.min_targets}
        meas_all = {tf: zscore_mean(Xm[:, [col_of[g] for g in ts if g in col_of]])[0]
                    for tf, ts in genesets.items()
                    if len([g for g in ts if g in col_of]) >= args.min_targets}
    elif args.score == "aucell":
        log("scoring with AUCell -- WARNING: decoupler's n_up defaults to ceil(0.05*nvar), which")
        log("  saturates the AUC to a constant for small gene sets. Verified against this")
        log("  dataset's own stored AUCell: agreement is good for large regulons (r 0.63-0.90)")
        log("  but 3/20 small regulons get EXACTLY zero variance. Do not trust small regulons here.")
        proj_all = aucell_scores(Xp, needed, genesets, args.min_targets)
        meas_all = aucell_scores(Xm, needed, genesets, args.min_targets)
    else:
        log("scoring with mean per-gene z (NOT invariant to per-spot scale -- see --score doc)")
        proj_all = {tf: zscore_mean(Xp[:, [col_of[g] for g in ts if g in col_of]])[0]
                    for tf, ts in genesets.items()
                    if len([g for g in ts if g in col_of]) >= args.min_targets}
        meas_all = {tf: zscore_mean(Xm[:, [col_of[g] for g in ts if g in col_of]])[0]
                    for tf, ts in genesets.items()
                    if len([g for g in ts if g in col_of]) >= args.min_targets}

    rows = []
    proj_scores_map, meas_scores_map = {}, {}
    for tf, targets in sorted(genesets.items()):
        idx = [col_of[g] for g in targets if g in col_of]
        if len(idx) < args.min_targets:
            continue
        pv, mv = proj_all.get(tf), meas_all.get(tf)
        if pv is None or mv is None:
            continue
        n_p = n_m = len(idx)
        rec = {
            "tf": tf,
            "n_targets_shared": len(idx),
            "n_targets_used_projected": n_p,
            "n_targets_used_measured": n_m,
            "morans_I_projected": morans_i(pv, W),
            "morans_I_measured": morans_i(mv, W),
            "spearman_proj_vs_meas": float(stats.spearmanr(pv, mv).statistic),
        }
        rec["delta_morans_I"] = rec["morans_I_projected"] - rec["morans_I_measured"]
        proj_scores_map[tf], meas_scores_map[tf] = pv, mv

        if chrom_df is not None:
            motifs = [m for m in motif_tf.get(tf.upper(), []) if m in chrom_df.columns]
            if not motifs:
                motifs = [c for c in chrom_df.columns if c.upper() == tf.upper()]
            if motifs:
                # If a TF has several motifs, pick the one with the strongest PROJECTED
                # agreement. That favours the projected layer, making the measured-layer
                # comparison conservative rather than flattering.
                best_m, best_r = None, -2.0
                for m in motifs:
                    r = stats.spearmanr(pv, chrom_df[m].to_numpy()).statistic
                    if np.isfinite(r) and r > best_r:
                        best_m, best_r = m, float(r)
                if best_m is not None:
                    mo = chrom_df[best_m].to_numpy()
                    rec.update({
                        "motif": best_m,
                        "n_motifs_for_tf": len(motifs),
                        "r_projected_vs_chromvar": best_r,
                        "r_measured_vs_chromvar": float(stats.spearmanr(mv, mo).statistic),
                    })
                    rec["drop"] = (rec["r_projected_vs_chromvar"]
                                   - rec["r_measured_vs_chromvar"])
        rows.append(rec)

    if not rows:
        raise SystemExit(f"no regulon had >= {args.min_targets} targets shared by both layers")
    df = pd.DataFrame(rows)

    # ---- report ----
    log("")
    log("=" * 78)
    log("MORAN'S I INFLATION FROM PROJECTION (per regulon, identical score, same graph)")
    log("=" * 78)
    log(f"  regulons evaluated        : {len(df)}")
    log(f"  median I, projected layer : {df['morans_I_projected'].median():+.4f}")
    log(f"  median I, measured layer  : {df['morans_I_measured'].median():+.4f}")
    log(f"  median delta (proj-meas)  : {df['delta_morans_I'].median():+.4f}")
    n_inf = int((df["delta_morans_I"] > 0).sum())
    log(f"  inflated by projection    : {n_inf}/{len(df)}")
    if len(df) >= 6:
        w = stats.wilcoxon(df["morans_I_projected"], df["morans_I_measured"])
        log(f"  Wilcoxon across regulons (n={len(df)}): stat={w.statistic:.1f}, p={w.pvalue:.3e}")
        log("    regulons are the replicate unit, NOT spots -- spot-wise tests would be")
        log("    pseudoreplication (paper_use_cases.md, statistical requirements)")
    log(f"  median Spearman(projected, measured) per regulon = "
        f"{df['spearman_proj_vs_meas'].median():+.3f}")

    # ---- degeneracy: do the regulon scores actually differ from one another? ----
    # Prompted by the raw table: morans_I_projected sits at ~0.85 for nearly every regulon
    # regardless of which genes it contains, while the measured layer spreads from +0.14 to
    # +0.89. If projected regulon scores are mutually near-identical, then that ~0.85 is a
    # property of the projection rather than of any regulon, and ranking regulons by spatial
    # structure on projected data is meaningless. Tangram's output is a weighted average of
    # snRNA profiles, so it can be near-low-rank -- in which case every gene set inherits the
    # same handful of spatial components.
    def degeneracy(score_map, tfs):
        S = np.column_stack([score_map[t] for t in tfs])
        R = np.abs(np.corrcoef(S, rowvar=False))
        off = R[np.triu_indices_from(R, k=1)]
        Sc = (S - S.mean(axis=0)) / np.where(S.std(axis=0) > 0, S.std(axis=0), 1)
        ev = np.linalg.svd(Sc, compute_uv=False) ** 2
        return {
            "median_abs_pairwise_r": float(np.median(off)),
            "frac_pairs_above_0.9": float((off > 0.9).mean()),
            "pc1_variance_explained": float(ev[0] / ev.sum()),
            "participation_ratio": float(ev.sum() ** 2 / (ev ** 2).sum()),
        }

    tfs = list(df["tf"])
    deg_p = degeneracy(proj_scores_map, tfs)
    deg_m = degeneracy(meas_scores_map, tfs)
    log("")
    log("=" * 78)
    log("DEGENERACY: are the regulon scores mutually distinguishable within each layer?")
    log("=" * 78)
    log(f"{'':34}{'projected':>12}{'measured':>12}")
    for k in ("median_abs_pairwise_r", "frac_pairs_above_0.9",
              "pc1_variance_explained", "participation_ratio"):
        log(f"  {k:32}{deg_p[k]:>12.3f}{deg_m[k]:>12.3f}")
    log("  participation_ratio is an effective count of independent score dimensions")
    log(f"  (out of {len(tfs)} regulons)")

    # size dependence of the inflation, which the raw table suggests is the whole story
    if len(df) >= 6:
        rho = stats.spearmanr(df["n_targets_shared"], df["delta_morans_I"])
        log("")
        log(f"  Spearman(regulon size, delta Moran's I) = {rho.statistic:+.3f} "
            f"(p={rho.pvalue:.2e})")
        small = df[df["n_targets_shared"] < 200]
        large = df[df["n_targets_shared"] >= 200]
        log(f"  median delta, regulons <200 targets (n={len(small)}): "
            f"{small['delta_morans_I'].median():+.4f}")
        log(f"  median delta, regulons >=200 targets (n={len(large)}): "
            f"{large['delta_morans_I'].median():+.4f}")

    arm = df.dropna(subset=["drop"]) if "drop" in df.columns else df.iloc[0:0]
    log("")
    log("=" * 78)
    log("CIRCULARITY CONTROL: does chromatin-expression agreement survive de-projection?")
    log("=" * 78)
    if arm.empty:
        log("  no TF had both a scorable regulon and a motif in the panel -- arm not evaluable")
    else:
        log(f"  TFs with a regulon AND a chromVAR motif: {len(arm)}")
        log(f"  median r(projected regulon, chromVAR) = {arm['r_projected_vs_chromvar'].median():+.3f}")
        log(f"  median r(measured  regulon, chromVAR) = {arm['r_measured_vs_chromvar'].median():+.3f}")
        log(f"  median drop                           = {arm['drop'].median():+.3f}")
        if len(arm) >= 6:
            w = stats.wilcoxon(arm["r_projected_vs_chromvar"], arm["r_measured_vs_chromvar"])
            log(f"  Wilcoxon across TFs (n={len(arm)}): stat={w.statistic:.1f}, p={w.pvalue:.3e}")
        log(f"  TFs retaining r>0.1 on the measured layer: "
            f"{int((arm['r_measured_vs_chromvar'] > 0.1).sum())}/{len(arm)}")
        med_m = arm["r_measured_vs_chromvar"].median()
        med_d = arm["drop"].median()
        log("")
        log("  INTERPRETATION")
        if med_m > 0.1 and med_d < 0.15:
            log("  -> agreement largely SURVIVES de-projection: the co-localisation is not")
            log("     merely a shared-mapping artefact, which is what PRIMARY 3 needs.")
        elif med_m <= 0.05:
            log("  -> agreement COLLAPSES once the shared projection is removed: the")
            log("     published-style claim would be artefact; the diagnostic IS the result.")
        else:
            log("  -> partial survival. Report the bound, not a binary verdict.")

    log("")
    show = [c for c in ("tf", "n_targets_shared", "morans_I_projected", "morans_I_measured",
                        "delta_morans_I", "spearman_proj_vs_meas", "motif",
                        "r_projected_vs_chromvar", "r_measured_vs_chromvar", "drop")
            if c in df.columns]
    with pd.option_context("display.width", 220, "display.max_columns", 30):
        print(df[show].to_string(index=False, float_format=lambda v: f"{v:+.4f}"), flush=True)

    payload = {
        "inputs": {"projected": projected, "measured": measured},
        "n_spots_shared": len(shared),
        "n_genes_universe": len(universe),
        "n_target_genes_read": len(needed),
        "gene_bridge": best,
        "estimator": args.score,
        "n_regulons_evaluated": int(len(df)),
        "morans": {
            "median_projected": float(df["morans_I_projected"].median()),
            "median_measured": float(df["morans_I_measured"].median()),
            "median_delta": float(df["delta_morans_I"].median()),
            "n_inflated": n_inf,
        },
        "degeneracy": {"projected": deg_p, "measured": deg_m},
        "chromatin_arm": ({} if arm.empty else {
            "n_tfs": int(len(arm)),
            "median_r_projected": float(arm["r_projected_vs_chromvar"].median()),
            "median_r_measured": float(arm["r_measured_vs_chromvar"].median()),
            "median_drop": float(arm["drop"].median()),
        }),
        "per_regulon": json.loads(df.to_json(orient="records")),
    }
    with open(args.out, "w") as fh:
        json.dump(payload, fh, indent=2)
    log(f"wrote {args.out}")
    return 0


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--inventory", action="store_true",
                    help="Report structure of both layers (metadata only) and exit.")
    ap.add_argument("--stage", action="store_true",
                    help=f"Copy inputs to {STAGE_DIR} (local overlay, not NFS) before reading.")
    ap.add_argument("--projected", default=PROJECTED,
                    help="Override the projected-layer path (e.g. a local staged copy).")
    ap.add_argument("--measured", default=MEASURED,
                    help="Override the measured-layer path. The copy under data/heart is an "
                         "unreadable NFS object; the intact original is the scratchpad's "
                         "heart_st_final.h5ad, byte-identical at 128,542,115 bytes.")
    ap.add_argument("--out", default="heart_projection_control_results.json")
    ap.add_argument("--score",
                    choices=("meanz-rownorm", "meanz-norm", "aucell", "meanz"),
                    default="meanz-rownorm",
                    help="Regulon estimator, applied identically to both layers. 'aucell' is "
                         "rank-based per spot and therefore invariant to the per-spot scale "
                         "factor the projected layer carries -- but decoupler's n_up default "
                         "saturates it for small gene sets, so prefer 'meanz-norm', which "
                         "normalises both layers to a common per-spot scale first. Plain "
                         "'meanz' is retained only to reproduce the artefact that skipping "
                         "normalisation produces.")
    ap.add_argument("--min-targets", type=int, default=10,
                    help="Skip regulons with fewer than this many targets shared by both layers.")
    args = ap.parse_args()

    projected, measured = args.projected, args.measured
    for p in (projected, measured):
        if not os.path.exists(p):
            raise SystemExit(f"not found: {p}")
    if args.stage:
        log(f"staging inputs to {STAGE_DIR} (one sequential NFS read each)")
        projected, measured = stage([projected, measured])

    if args.inventory:
        inventory(projected, measured)
        return 0
    return control(args, projected, measured)


if __name__ == "__main__":
    sys.exit(main())
