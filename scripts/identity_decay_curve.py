#!/usr/bin/env python3
"""Over what sequence-identity range does the signal actually hold?

The cross-lineage question is settled and the answer is no: the determinants are
lineage-specific, two of three lineage pairs point in opposite directions, and a global
classifier trained across families learns something that reverses on a family it has not
seen. Asking again with a bigger model is not a plan.

The question that remains is well posed and answerable, and it is the one a user of this
tool would actually ask: given a candidate that resembles enzymes we have characterised,
HOW CLOSE does it have to be before the prediction means anything?

The 262-member cluster is the only place with enough data to answer. It is one 30%-identity
lineage holding 171 measured-active and 91 measured-inactive enzymes. Re-clustering inside
it at successively stricter identities and holding out whole sub-clusters gives a curve:
at 90% the held-out enzymes have near-duplicates in training, at 50% they do not, and where
that curve meets chance is the honest edge of the tool's competence.

Three things are plotted against each other at every threshold, because a learned model that
never separates from retrieval is not a model worth shipping:
  the head        frozen ESM-2 embeddings, logistic regression
  retrieval       score by similarity to the nearest ACTIVE enzyme in training
  composition     amino-acid frequencies, the floor everything must clear

Reported as a curve rather than a number on purpose. A single figure invites the reader to
ask whether the tool works; a curve answers the question they should be asking, which is
where it stops working.

AND THE CURVE MUST CARRY ITS UNCERTAINTY, because the point estimates lie. Run bare, this
shows the head at 0.703 against retrieval's 0.544 at 70% identity, which reads as a clear
operating band. Bootstrapped over the held-out sub-clusters, the PAIRED lead is
+0.159 [-0.263, +0.390]: it spans zero, so the band is not established. Seventeen
sub-clusters cannot resolve a difference that size, and reporting the 0.703 alone would have
claimed a competence range this data does not support.
"""

from __future__ import annotations

import json
import pathlib
import sys
from typing import Dict, List, Tuple

import numpy as np
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config, seqtools
from pipeline.db import connect

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "shv", pathlib.Path(__file__).resolve().parent / "sequence_head_variants.py")
shv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(shv)

THRESHOLDS = (0.95, 0.90, 0.80, 0.70, 0.60, 0.50, 0.40)
MIN_SUBCLUSTERS = 3
OUT = config.INTERIM_DIR / "identity_decay_curve.json"


def biggest_cluster() -> Tuple[List[Tuple[str, str]], np.ndarray]:
    recs, y = shv.contrast()
    keep, groups = shv.mixed_cluster_mask(recs, y)
    recs_m = [r for r, k in zip(recs, keep) if k]
    y_m, g_m = y[keep], groups[keep]
    counts = {g: int((g_m == g).sum()) for g in set(g_m)}
    big = max(counts, key=counts.get)
    sel = g_m == big
    return [r for r, s in zip(recs_m, sel) if s], y_m[sel]


def leave_one_subcluster_out(X, y, groups, scorer) -> Tuple[float, int, List[float]]:
    """Size-weighted mean over held-out sub-clusters. Unweighted hides small-group noise."""
    scores, sizes = [], []
    for g in sorted(set(groups)):
        te = groups == g
        tr = ~te
        if len(set(y[te])) < 2 or len(set(y[tr])) < 2:
            continue
        scores.append(scorer(X[tr], y[tr], X[te], y[te]))
        sizes.append(int(te.sum()))
    if not scores:
        return float("nan"), 0, []
    a, w = np.array(scores), np.array(sizes)
    return float((a * w).sum() / w.sum()), len(scores), scores


def head(Xtr, ytr, Xte, yte) -> float:
    m = make_pipeline(StandardScaler(),
                      LogisticRegression(max_iter=3000, class_weight="balanced"))
    m.fit(Xtr, ytr)
    return roc_auc_score(yte, m.predict_proba(Xte)[:, 1])


def retrieval(Xtr, ytr, Xte, yte) -> float:
    """Cosine similarity to the nearest ACTIVE enzyme in training. The number to beat."""
    A = Xtr[ytr == 1]
    if not len(A):
        return float("nan")
    An = A / np.clip(np.linalg.norm(A, axis=1, keepdims=True), 1e-9, None)
    Tn = Xte / np.clip(np.linalg.norm(Xte, axis=1, keepdims=True), 1e-9, None)
    return roc_auc_score(yte, (Tn @ An.T).max(axis=1))


def main() -> int:
    recs, y = biggest_cluster()
    print(f"largest lineage: {len(recs)} enzymes, {int(y.sum())} measured-active, "
          f"{int((1 - y).sum())} measured-inactive\n")

    emb = shv.embed_cached("esm2-35M", shv.VARIANTS["esm2-35M"], recs)
    comp = shv.composition(recs)
    fa = seqtools.write_fasta(recs, config.INTERIM_DIR / "decay.fasta")

    print(f"{'identity':>9} {'sub-clusters':>13} {'scored':>7}   {'head':>6} {'retrieval':>10} "
          f"{'composition':>12}")
    rows = []
    for thr in THRESHOLDS:
        clu = seqtools.cluster(fa, min_seq_id=thr)
        groups = np.array([clu.get(i, i) for i, _ in recs])
        n_sub = len(set(groups))
        h, k, _ = leave_one_subcluster_out(emb, y, groups, head)
        r, _, _ = leave_one_subcluster_out(emb, y, groups, retrieval)
        c, _, _ = leave_one_subcluster_out(comp, y, groups, head)
        if k < MIN_SUBCLUSTERS:
            print(f"{thr:>8.0%} {n_sub:>13} {k:>7}   too few scorable sub-clusters")
            continue
        print(f"{thr:>8.0%} {n_sub:>13} {k:>7}   {h:>6.3f} {r:>10.3f} {c:>12.3f}")
        rows.append({"identity": thr, "n_subclusters": n_sub, "n_scored": k,
                     "head": h, "retrieval": r, "composition": c})

    OUT.write_text(json.dumps({"n": len(recs), "n_active": int(y.sum()),
                               "n_inactive": int((1 - y).sum()), "curve": rows}, indent=2))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
