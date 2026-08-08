#!/usr/bin/env python3
"""Try to improve the sequence head, measured against one fixed contrast.

The sequence head is the only signal in this project that survived measured labels:
0.593 +/- 0.078 on the within-family contrast restricted to mixed clusters, against an
amino-acid composition baseline of 0.440. Geometry, on the same labels, is at 0.398.

So this exists to answer one question -- can that 0.593 be improved -- and to answer it in
a way that cannot flatter itself. Every variant is scored on exactly the same enzymes, the
same 30% identity clusters, the same folds and the same seeds, with only the REPRESENTATION
changing. Composition is recomputed alongside each one, because a representation that lifts
the head and the baseline equally has added nothing.

The contrast is deliberately the hardest available: measured PET-active against measured
PET-inactive, both polyesterases, restricted to the clusters that contain both classes so
the model cannot win by recognising a lineage. That restriction is why the numbers here are
so much lower than the 0.98 out-of-family figure, and it is the only one worth optimising.

Variants:
  esm2-35M    the current representation, 480 dimensions. The baseline.
  esm2-150M   640 dimensions.
  esm2-650M   1280 dimensions. Roughly eighteen times the parameters of the current model.

A larger protein language model is the obvious first thing to try and the least
interesting: if the answer is that within-family PET activity is not written in the
sequence in a way a frozen embedding exposes, no amount of extra dimensions will find it,
and that is worth establishing before anything cleverer is attempted.
"""

from __future__ import annotations

import json
import pathlib
import sys
import time
from typing import Dict, List, Tuple

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config, seqtools
from pipeline.db import connect
from pipeline.embed import esm

VARIANTS = {
    "esm2-35M": "facebook/esm2_t12_35M_UR50D",
    "esm2-150M": "facebook/esm2_t30_150M_UR50D",
    "esm2-650M": "facebook/esm2_t33_650M_UR50D",
}
SEEDS = 10
OUT = config.INTERIM_DIR / "sequence_head_variants.json"
CACHE = config.INTERIM_DIR / "variant_embeddings"

AA = "ACDEFGHIKLMNPQRSTVWY"


def contrast() -> Tuple[List[Tuple[str, str]], np.ndarray]:
    """Measured PET-active against measured PET-inactive polyesterases."""
    marks = ",".join("?" * len(config.MEASURED_TIERS))
    with connect() as c:
        pos = c.execute(
            f"SELECT enzyme_id, sequence FROM characterised_enzymes "
            f"WHERE is_positive=1 AND source_ref IN ({marks}) AND excluded_from_training=0 "
            f"  AND sequence IS NOT NULL AND family!='mhetase_like'",
            config.MEASURED_TIERS).fetchall()
        neg = c.execute(
            "SELECT enzyme_id, sequence FROM characterised_enzymes "
            "WHERE within_family_basis='measured-inactive' AND excluded_from_training=0 "
            "  AND sequence IS NOT NULL").fetchall()
    recs = [tuple(r) for r in pos] + [tuple(r) for r in neg]
    y = np.array([1] * len(pos) + [0] * len(neg))
    return recs, y


def mixed_cluster_mask(recs, y) -> Tuple[np.ndarray, np.ndarray]:
    """Keep only clusters holding BOTH classes, and return the grouping."""
    fa = seqtools.write_fasta(recs, config.INTERIM_DIR / "variants.fasta")
    clu = seqtools.cluster(fa, min_seq_id=0.30)
    groups = np.array([clu.get(i, i) for i, _ in recs])
    shared = set(groups[y == 1]) & set(groups[y == 0])
    keep = np.array([g in shared for g in groups])
    return keep, groups


def grouped_auc(X: np.ndarray, y: np.ndarray, groups: np.ndarray) -> Tuple[float, float]:
    """Leave-one-cluster-out, reported as the mean and spread ACROSS clusters.

    Was k-fold over clusters with a standard deviation taken across random seeds, which
    was measuring nothing. The clusters here are sized 262, 57, 38, 6, 5, 4, 4, 3, 3, 2 --
    one holds 68% of the data -- so with four folds the partition is nearly forced: ten
    seeds produced two distinct splits, and the "+/- 0.000" that came out of it was
    reporting that fact rather than a stable estimate.

    Holding out one cluster at a time removes the dependence on a partition entirely, and
    the spread is then across lineages, which is the quantity that actually matters here:
    it says whether the signal transfers to a family the model has not seen, and a mean
    that hides one cluster carrying everything is worth distrusting.
    """
    scores, sizes = [], []
    for g in sorted(set(groups)):
        te = groups == g
        tr = ~te
        if len(set(y[te])) < 2 or len(set(y[tr])) < 2:
            continue
        model = make_pipeline(StandardScaler(),
                              LogisticRegression(max_iter=3000, class_weight="balanced"))
        model.fit(X[tr], y[tr])
        scores.append(roc_auc_score(y[te], model.predict_proba(X[te])[:, 1]))
        sizes.append(int(te.sum()))
    if not scores:
        return float("nan"), float("nan")
    a, w = np.array(scores), np.array(sizes)
    # The UNWEIGHTED mean is the trap. Seven of these ten clusters hold two to six enzymes,
    # where a leave-one-out AUC can only come out near 0 or near 1, and averaging those
    # seven coin flips alongside the three clusters large enough to measure anything is
    # what produced a headline of 0.625 from three real values of 0.443, 0.429 and 0.574.
    # Both are returned; the size-weighted one is the honest summary.
    return float((a * w).sum() / w.sum()), float(a.std())


def composition(recs) -> np.ndarray:
    return np.array([[s.count(x) / max(len(s), 1) for x in AA] for _, s in recs])


def embed_cached(name: str, model_id: str, recs) -> np.ndarray:
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{name}.npz"
    if path.exists():
        d = np.load(path, allow_pickle=True)
        by_id = {i: v for i, v in zip(d["ids"], d["X"])}
        if all(i in by_id for i, _ in recs):
            print(f"    (cached)")
            return np.vstack([by_id[i] for i, _ in recs])
    t0 = time.time()
    ids, X, rep = esm.embed(recs, batch_size=4, model_name=model_id, progress_every=100)
    print(f"    embedded {rep['n']} at dim {rep['dim']} in {time.time()-t0:.0f}s")
    np.savez_compressed(path, ids=np.array(ids, dtype=object), X=X)
    by_id = {i: v for i, v in zip(ids, X)}
    return np.vstack([by_id[i] for i, _ in recs])


def main() -> int:
    recs, y = contrast()
    keep, groups = mixed_cluster_mask(recs, y)
    recs_m = [r for r, k in zip(recs, keep) if k]
    y_m, g_m = y[keep], groups[keep]
    print(f"full contrast : {len(recs)} enzymes, {int(y.sum())} active, {int((1-y).sum())} inactive")
    print(f"mixed clusters: {len(recs_m)} enzymes, {int(y_m.sum())} active, "
          f"{int((1-y_m).sum())} inactive, {len(set(g_m))} clusters\n")

    comp_mean, comp_sd = grouped_auc(composition(recs_m), y_m, g_m)
    print(f"amino-acid composition baseline: {comp_mean:.3f} +/- {comp_sd:.3f}\n")

    results = {"composition": {"auc": comp_mean, "sd": comp_sd},
               "n": len(recs_m), "n_clusters": len(set(g_m)), "seeds": SEEDS}
    for name, model_id in VARIANTS.items():
        print(f"{name}:")
        try:
            X = embed_cached(name, model_id, recs)
        except Exception as exc:
            print(f"    FAILED: {exc}")
            results[name] = {"error": str(exc)}
            continue
        Xm = X[keep]
        mean, sd = grouped_auc(Xm, y_m, g_m)
        results[name] = {"auc_size_weighted": mean, "sd_across_clusters": sd,
                         "dim": int(X.shape[1]), "over_composition": mean - comp_mean}
        print(f"    size-weighted AUC {mean:.3f}   (spread across clusters {sd:.3f})   "
              f"dim {X.shape[1]}   over composition {mean - comp_mean:+.3f}\n")

    OUT.write_text(json.dumps(results, indent=2))
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
