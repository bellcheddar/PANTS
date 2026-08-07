#!/usr/bin/env python3
"""Does active-site geometry track PET activity, when the labels were MEASURED?

Every earlier version of this test used negatives that mean "not reported active on PET",
which cannot separate an enzyme that was assayed and failed from one nobody assayed. On
those labels the raw signal looked convincing -- cleft depth at AUC 0.808, p 1.7e-07 -- and
collapsed to 0.534 under cluster-grouped splitting, which was read as an underpowered
evaluation.

The screen ingests supply the label that was missing: proteins expressed, assayed under a
published protocol, and found to release no product. This asks the same question of the
same features against those.

Restricted to ESMFold structures on BOTH sides. Mixing coordinate sources is the confound
established by the paired comparison in structure_source_confound.py -- the oxyanion hole
differs systematically between a crystal structure and a model of the same protein, n=51,
p 1.1e-06 -- so an all-predicted contrast is the only internally consistent one available,
and it is stated here rather than left for a reader to notice.
"""

from __future__ import annotations

import json
import pathlib
import sys

import numpy as np
from scipy import stats
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedGroupKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config, seqtools
from pipeline.db import connect

FEATURES = ["ser_og_his_ne2_dist_A", "his_nd1_asp_od_dist_A", "ser_his_asp_angle_deg",
            "oxyanion_n1_dist_A", "oxyanion_n2_dist_A", "oxyanion_n2_angle_deg",
            "cleft_width_A", "cleft_depth_A", "n_cleft_residues"]

OUT = config.INTERIM_DIR / "geometry_measured_labels.json"
SEEDS = 10


def load() -> list:
    with connect() as c:
        rows = [dict(r) for r in c.execute(f"""
            SELECT ce.enzyme_id, ce.sequence, ce.within_family_basis,
                   {', '.join('rg.' + f for f in FEATURES)}
            FROM reference_geometry rg
            JOIN reference_structures rs USING(enzyme_id)
            JOIN characterised_enzymes ce USING(enzyme_id)
            WHERE rs.source = 'esmfold'
              AND rg.triad_ser_resnum IS NOT NULL AND rg.cleft_depth_A IS NOT NULL
              AND (ce.is_positive = 1 OR ce.within_family_basis = 'measured-inactive')""")]
    return [r for r in rows if all(r[f] is not None for f in FEATURES)]


def main() -> int:
    rows = load()
    y = np.array([0 if r["within_family_basis"] == "measured-inactive" else 1 for r in rows])
    X = np.array([[r[f] for f in FEATURES] for r in rows], dtype=float)
    fa = seqtools.write_fasta([(r["enzyme_id"], r["sequence"]) for r in rows],
                              config.INTERIM_DIR / "geom_measured.fasta")
    clu = seqtools.cluster(fa, min_seq_id=0.3)
    groups = np.array([clu.get(r["enzyme_id"], r["enzyme_id"]) for r in rows])

    shared = set(groups[y == 0]) & set(groups[y == 1])
    print(f"{len(rows)} ESMFold structures: {int(y.sum())} PET-active, "
          f"{int((1 - y).sum())} measured-inactive")
    print(f"{len(set(groups))} clusters, {len(shared)} holding both classes\n")

    per_feature = {}
    print("raw feature differences (no splitting -- the optimistic reading)")
    for i, f in enumerate(FEATURES):
        auc = roc_auc_score(y, X[:, i])
        p = float(stats.mannwhitneyu(X[y == 1, i], X[y == 0, i]).pvalue)
        per_feature[f] = {"auc": auc, "p": p}
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""
        print(f"   {f:24s} AUC {auc:.3f}  p {p:.2g} {star}")

    means = []
    for seed in range(SEEDS):
        cv = StratifiedGroupKFold(n_splits=5, shuffle=True, random_state=seed)
        model = make_pipeline(StandardScaler(),
                              LogisticRegression(max_iter=2000, class_weight="balanced"))
        scores = []
        for tr, te in cv.split(X, y, groups):
            if len(set(y[tr])) < 2 or len(set(y[te])) < 2:
                continue
            model.fit(X[tr], y[tr])
            scores.append(roc_auc_score(y[te], model.predict_proba(X[te])[:, 1]))
        if scores:
            means.append(float(np.mean(scores)))
    means_a = np.array(means)
    print(f"\ncluster-grouped, all features, over {SEEDS} seeds")
    print(f"   AUC {means_a.mean():.3f}   range {means_a.min():.3f}-{means_a.max():.3f}")
    print("\n   on INFERRED labels the same test gave 0.534 +/- 0.173,")
    print("   and cleft depth alone gave 0.808 at p 1.7e-07.")

    OUT.write_text(json.dumps({
        "n": len(rows), "n_positive": int(y.sum()), "n_negative": int((1 - y).sum()),
        "n_clusters": len(set(groups)), "n_shared_clusters": len(shared),
        "coordinate_source": "esmfold only, both classes",
        "grouped_auc_mean": float(means_a.mean()),
        "grouped_auc_min": float(means_a.min()), "grouped_auc_max": float(means_a.max()),
        "seeds": SEEDS, "per_feature": per_feature,
    }, indent=2))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
