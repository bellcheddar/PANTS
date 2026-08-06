#!/usr/bin/env python3
"""Re-ask whether active-site geometry tracks PET activity, on the finished structure set.

The finding this supersedes was measured on 131 AlphaFold models: raw feature differences
looked convincing and physically sensible (cleft depth AUC 0.819, p<0.001), and under
cluster-grouped splitting it collapsed to 0.533 +/- 0.185, indistinguishable from the
sequence head's 0.493. The conclusion drawn was that the signature is largely cluster
structure and the evaluation is label-limited rather than method-limited.

That conclusion deserves re-testing rather than re-asserting, for two reasons. The set is
now 344 rather than 131, and a bigger benchmark has already overturned a "flat" reading
once on a sibling project. And 57 of these are crystal structures rather than models,
which matters more than the count: the previous run measured geometry off predictions in
exactly the flexible loops that gate the active site, which is where a prediction is least
reliable and where this project looks.

So three questions, not one:

  1. does the raw signal reproduce at 2.6x the size?
  2. does it survive cluster-grouped splitting this time?
  3. do experimental structures give a different answer from predicted ones? If the
     signature is real, the crystal structures should show it at least as strongly. If it
     is a modelling artefact, it should weaken where the coordinates are measured.

Reads reference_geometry rather than re-downloading and re-measuring: those rows were
written by the same geometry.measure() the old script called, so this compares like with
like while adding the experimental structures the old one could not see.
"""

from __future__ import annotations

import json
import pathlib
import sys
from typing import Dict, List, Optional

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config, seqtools
from pipeline.db import connect

from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

FEATURES = ["ser_og_his_ne2_dist_A", "his_nd1_asp_od_dist_A", "ser_his_asp_angle_deg",
            "oxyanion_n1_dist_A", "oxyanion_n2_dist_A", "oxyanion_n2_angle_deg",
            "cleft_width_A", "cleft_depth_A", "n_cleft_residues"]

# Same contrast the embeddings failed: PET-active polyesterases against polyesterases
# measured on a different plastic. Not "PETase versus not-a-hydrolase", which is easy and
# already answered at AUC 0.975.
SQL = f"""
SELECT e.enzyme_id, e.sequence, r.source, r.plddt_mean, r.resolution_A,
       e.is_positive, e.within_family_basis,
       {', '.join('g.' + f for f in FEATURES)}
FROM reference_geometry g
JOIN reference_structures r USING(enzyme_id)
JOIN characterised_enzymes e USING(enzyme_id)
WHERE g.triad_ser_resnum IS NOT NULL AND g.cleft_depth_A IS NOT NULL
  AND (e.is_positive = 1 OR e.within_family_basis IS NOT NULL)
"""


def load() -> List[dict]:
    with connect() as c:
        rows = [dict(r) for r in c.execute(SQL)]
    for r in rows:
        # A within-family negative is a negative even when the source row also carries
        # is_positive, because the basis column is the measured statement about PET.
        r["label"] = 0 if r["within_family_basis"] else 1
    return [r for r in rows if all(r[f] is not None for f in FEATURES)]


def clusters(rows: List[dict]) -> Dict[str, str]:
    """30% identity clusters, so a fold is never split across train and test.

    Same call and the same threshold eval_within_family.py uses, so the two evaluations
    group on identical units and their numbers are comparable.
    """
    recs = [(r["enzyme_id"], r["sequence"]) for r in rows if r["sequence"]]
    fa = seqtools.write_fasta(recs, config.INTERIM_DIR / "geom_v2.fasta")
    return seqtools.cluster(fa, min_seq_id=0.3)


def auc_ci(y: np.ndarray, x: np.ndarray) -> tuple:
    """AUC with a Mann-Whitney p, direction-free: report |AUC| against a coin flip."""
    if len(set(y)) < 2:
        return None, None
    a = roc_auc_score(y, x)
    u = stats.mannwhitneyu(x[y == 1], x[y == 0], alternative="two-sided")
    return a, u.pvalue


def grouped_auc(X: np.ndarray, y: np.ndarray, groups: np.ndarray,
                n_splits: int = 5) -> Optional[tuple]:
    """Cluster-grouped cross-validated AUC, per fold as well as pooled.

    Per-fold values are printed because the pooled mean hid the variance last time: a
    0.533 +/- 0.185 mean is a different statement from five folds that all sat near 0.53.
    """
    n_neg_groups = len(set(groups[y == 0]))
    if n_neg_groups < n_splits:
        n_splits = max(2, n_neg_groups)
    if n_neg_groups < 2:
        return None
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=0)
    model = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000,
                                                               class_weight="balanced"))
    scores = []
    for tr, te in cv.split(X, y, groups):
        if len(set(y[tr])) < 2 or len(set(y[te])) < 2:
            continue
        model.fit(X[tr], y[tr])
        scores.append(roc_auc_score(y[te], model.predict_proba(X[te])[:, 1]))
    if not scores:
        return None
    return float(np.mean(scores)), float(np.std(scores)), scores, n_splits


def report(name: str, rows: List[dict], cl: Dict[str, int]) -> None:
    y = np.array([r["label"] for r in rows])
    n_pos, n_neg = int((y == 1).sum()), int((y == 0).sum())
    # An enzyme MMseqs2 did not place is its own cluster, never a shared bucket: pooling
    # the unplaced would invent a group that leaks across the split.
    groups = np.array([cl.get(r["enzyme_id"], r["enzyme_id"]) for r in rows])
    print(f"\n=== {name} ===")
    print(f"{len(rows)} enzymes: {n_pos} PET-active, {n_neg} within-family negatives")
    print(f"{len(set(groups))} clusters, {len(set(groups[y == 0]))} of them holding a negative")
    if n_neg < 2:
        print("  not evaluable: fewer than two negatives")
        return

    print("\n  raw feature differences (no splitting -- the optimistic reading)")
    for f in FEATURES:
        x = np.array([r[f] for r in rows], dtype=float)
        a, p = auc_ci(y, x)
        if a is None:
            continue
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        print(f"    {f:24s} AUC {a:.3f}  p {p:.2g} {star}")

    X = np.array([[r[f] for f in FEATURES] for r in rows], dtype=float)
    got = grouped_auc(X, y, groups)
    print("\n  cluster-grouped, all features (the honest reading)")
    if got is None:
        print("    not evaluable: too few independent clusters hold a negative")
    else:
        mean, sd, scores, k = got
        print(f"    AUC {mean:.3f} +/- {sd:.3f} over {k} folds")
        print(f"    per fold: {', '.join(f'{s:.3f}' for s in scores)}")

    plddt = [r["plddt_mean"] for r in rows if r["plddt_mean"] is not None]
    if plddt:
        pos = [r["plddt_mean"] for r in rows if r["label"] == 1 and r["plddt_mean"]]
        neg = [r["plddt_mean"] for r in rows if r["label"] == 0 and r["plddt_mean"]]
        if pos and neg:
            # A geometric difference that tracks model confidence is an artefact, so the
            # confidence gap is reported next to the result rather than checked once.
            print(f"\n  mean pLDDT: active {np.mean(pos):.1f}, negative {np.mean(neg):.1f} "
                  f"(a gap here would make the result a modelling artefact)")


def main() -> int:
    rows = load()
    if not rows:
        print("no rows with complete geometry")
        return 1
    cl = clusters(rows)

    report("all structures", rows, cl)
    report("experimental only (crystal structures)",
           [r for r in rows if r["source"] == "pdb"], cl)
    report("predicted only (AlphaFold and ESMFold)",
           [r for r in rows if r["source"] != "pdb"], cl)

    out = config.INTERIM_DIR / "geometry_vs_activity_v2.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps({"n": len(rows), "features": FEATURES,
                               "by_source": {s: sum(1 for r in rows if r["source"] == s)
                                             for s in {r["source"] for r in rows}}}, indent=2))
    print(f"\nwrote {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
