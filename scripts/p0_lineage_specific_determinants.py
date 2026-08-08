#!/usr/bin/env python3
"""Go/no-go: do different lineages agree about what makes a PETase?

Everything else waits on this. A global sequence-to-activity model can only exist if the
features that separate active from inactive point in roughly the same direction inside one
lineage as inside another. If the directions are unrelated, the determinants are
lineage-specific and no objective, architecture or amount of regularisation will produce a
transferable classifier -- and that is a finding about the biology rather than a failure of
the modelling.

The design has no train/test split anywhere, so it cannot leak. Fit the same model
separately inside each large cluster; compare the fitted coefficient vectors between
clusters.

THE COMPARISON NEEDS ITS OWN NULL, and getting that right is the whole test. A cosine
similarity of 0.3 between two coefficient vectors means nothing on its own: with 9 features
and 57 enzymes the vectors are noisy, so even two samples from an IDENTICAL underlying
process would not agree perfectly. So the reference is not zero. It is the agreement two
bootstrap replicates of the SAME cluster show with each other -- an empirical ceiling on
what agreement can look like when the true direction is by construction the same. Between-
cluster agreement is then read against that ceiling, not against zero.

Read it as:
  between-cluster cosine near the within-cluster ceiling  -> lineages agree, global model alive
  between-cluster cosine near zero, ceiling high          -> lineage-specific determinants
  ceiling itself low                                      -> too noisy to answer either way,
                                                             and no other proposal here is
                                                             answerable either
"""

from __future__ import annotations

import json
import pathlib
import sys
from typing import Dict, List, Tuple

import numpy as np
from sklearn.decomposition import PCA
from sklearn.linear_model import LogisticRegression
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config
from pipeline.db import connect

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "shv", pathlib.Path(__file__).resolve().parent / "sequence_head_variants.py")
shv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(shv)

GEOM = ["ser_og_his_ne2_dist_A", "his_nd1_asp_od_dist_A", "ser_his_asp_angle_deg",
        "oxyanion_n1_dist_A", "oxyanion_n2_dist_A", "oxyanion_n2_angle_deg",
        "cleft_width_A", "cleft_depth_A", "n_cleft_residues"]
BOOTSTRAPS = 200
MIN_CLUSTER = 30
OUT = config.INTERIM_DIR / "p0_lineage_determinants.json"


def fit_direction(X: np.ndarray, y: np.ndarray, seed: int = 0) -> np.ndarray:
    """Unit-length coefficient vector, features standardised so they are comparable.

    Standardisation is per fit and therefore per cluster, which is the point: the question
    is about DIRECTION within a lineage, not about absolute scales that differ between them.
    """
    model = make_pipeline(StandardScaler(),
                          LogisticRegression(max_iter=5000, class_weight="balanced", C=1.0,
                                             random_state=seed))
    model.fit(X, y)
    w = model[-1].coef_.ravel()
    n = np.linalg.norm(w)
    return w / n if n else w


def bootstrap_ceiling(X: np.ndarray, y: np.ndarray, rng: np.random.Generator) -> List[float]:
    """Cosine between two independent bootstrap resamples of the SAME cluster.

    This is the empirical ceiling: how much two fits agree when the underlying direction is
    identical and only sampling noise differs.
    """
    out = []
    n = len(y)
    for _ in range(BOOTSTRAPS):
        cos = []
        for _ in range(2):
            for _try in range(20):
                idx = rng.integers(0, n, n)
                if len(set(y[idx])) == 2:
                    break
            else:
                break
            cos.append(fit_direction(X[idx], y[idx]))
        if len(cos) == 2:
            out.append(float(cos[0] @ cos[1]))
    return out


def between(a: Tuple[np.ndarray, np.ndarray], b: Tuple[np.ndarray, np.ndarray],
            rng: np.random.Generator) -> List[float]:
    """Cosine between bootstrap fits of two DIFFERENT clusters."""
    Xa, ya = a
    Xb, yb = b
    out = []
    for _ in range(BOOTSTRAPS):
        vs = []
        for X, y in ((Xa, ya), (Xb, yb)):
            n = len(y)
            for _try in range(20):
                idx = rng.integers(0, n, n)
                if len(set(y[idx])) == 2:
                    break
            else:
                break
            vs.append(fit_direction(X[idx], y[idx]))
        if len(vs) == 2:
            out.append(float(vs[0] @ vs[1]))
    return out


def load_features() -> Tuple[Dict[str, np.ndarray], Dict[str, np.ndarray], np.ndarray, list]:
    recs, y = shv.contrast()
    keep, groups = shv.mixed_cluster_mask(recs, y)
    recs_m = [r for r, k in zip(recs, keep) if k]
    y_m, g_m = y[keep], groups[keep]
    ids = [e for e, _ in recs_m]

    with connect() as c:
        geo_rows = {r[0]: r[1:] for r in c.execute(
            "SELECT rg.enzyme_id, " + ", ".join("rg." + f for f in GEOM) +
            " FROM reference_geometry rg JOIN reference_structures rs USING(enzyme_id) "
            " WHERE rs.source='esmfold' AND rg.triad_ser_resnum IS NOT NULL")}

    emb = shv.embed_cached("esm2-35M", shv.VARIANTS["esm2-35M"], recs_m)
    return geo_rows, {e: v for e, v in zip(ids, emb)}, y_m, list(zip(ids, g_m))


def main() -> int:
    geo_rows, emb, y_m, id_grp = load_features()
    rng = np.random.default_rng(0)
    ids = [i for i, _ in id_grp]
    grp = np.array([g for _, g in id_grp])
    label = {i: int(v) for i, v in zip(ids, y_m)}

    big = [g for g in sorted(set(grp)) if (grp == g).sum() >= MIN_CLUSTER]
    print(f"clusters with at least {MIN_CLUSTER} enzymes: {len(big)}")

    report = {}
    for name, getter, dim in (
        ("geometry", lambda e: geo_rows.get(e), len(GEOM)),
        ("esm2-35M (PCA 20)", lambda e: emb.get(e), 20),
    ):
        print(f"\n=== {name} ===")
        blocks = {}
        for g in big:
            members = [i for i, gg in id_grp if gg == g and getter(i) is not None
                       and not any(v is None for v in np.atleast_1d(getter(i)))]
            ys = np.array([label[i] for i in members])
            if len(set(ys)) < 2 or len(members) < MIN_CLUSTER:
                print(f"  {g[:12]:<14} {len(members):>4} usable — skipped (single class or too few)")
                continue
            X = np.array([np.asarray(getter(i), dtype=float) for i in members])
            blocks[g] = (X, ys)
            print(f"  {g[:12]:<14} {len(members):>4} enzymes, {int(ys.sum())} active, "
                  f"{int((1-ys).sum())} inactive")
        if len(blocks) < 2:
            print("  fewer than two usable clusters — not answerable")
            continue

        if dim < next(iter(blocks.values()))[0].shape[1]:
            allX = np.vstack([X for X, _ in blocks.values()])
            pca = PCA(n_components=dim, random_state=0).fit(allX)
            blocks = {g: (pca.transform(X), y) for g, (X, y) in blocks.items()}

        print(f"\n  within-cluster bootstrap ceiling (identical true direction):")
        ceilings = {}
        for g, (X, y) in blocks.items():
            c = bootstrap_ceiling(X, y, rng)
            ceilings[g] = c
            print(f"    {g[:12]:<14} cosine {np.mean(c):+.3f}  "
                  f"[{np.percentile(c,5):+.3f}, {np.percentile(c,95):+.3f}]")

        print(f"\n  between-cluster agreement:")
        keys = list(blocks)
        pairs = {}
        for i in range(len(keys)):
            for j in range(i + 1, len(keys)):
                b = between(blocks[keys[i]], blocks[keys[j]], rng)
                pairs[f"{keys[i][:10]} vs {keys[j][:10]}"] = b
                ceil = 0.5 * (np.mean(ceilings[keys[i]]) + np.mean(ceilings[keys[j]]))
                print(f"    {keys[i][:10]} vs {keys[j][:10]}: cosine {np.mean(b):+.3f}  "
                      f"[{np.percentile(b,5):+.3f}, {np.percentile(b,95):+.3f}]   "
                      f"ceiling {ceil:+.3f}")
        report[name] = {
            "ceiling": {g: float(np.mean(v)) for g, v in ceilings.items()},
            "between": {k: {"mean": float(np.mean(v)),
                            "p05": float(np.percentile(v, 5)),
                            "p95": float(np.percentile(v, 95))} for k, v in pairs.items()},
        }

    OUT.write_text(json.dumps(report, indent=2))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
