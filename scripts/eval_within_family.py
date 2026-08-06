"""Does the head discriminate PET activity, or does it just recognise the fold?

This is the question limitation 8 records and the one that decides whether PANTS is more
than a well-built homology search. The existing headline (AUC 0.976) is measured against
hard negatives drawn from OTHER alpha/beta-hydrolase families, and the near-miss test
scored 1.000 against a set that is entirely one ESTHER `Cutinase` family. Both ask "is
this a polyesterase?", which the annotation pipeline already answers.

The interesting question is the next one in: **among enzymes that are already
polyesterases, which act on PET?**

So the same head is evaluated against three negative regimes of increasing difficulty,
everything else held fixed:

  out-of-family   hard negatives from other alpha/beta-hydrolase families
  near-miss       ESTHER Cutinase-family members, no measured PET activity
  within-family   PAZy enzymes MEASURED on another plastic (PA, PUR, PLA, PBAT, PHA),
                  restricted to those sharing a 30% cluster with a PET-active enzyme

and a fourth, strictest variant: within-family restricted to the mixed clusters alone, so
every negative has a PET-active enzyme inside its own cluster and separating them cannot
be done by recognising which family the sequence belongs to.

A composition-only baseline runs alongside each. That baseline is not decoration: it is
what a result has to beat to mean anything, and on this project it has already caught one
apparently strong number that was really amino-acid composition.

**The caveat that travels with the within-family numbers.** PAZy records only positive
substrate associations, so "PET not listed" conflates *inactive* with *never assayed*.
Some fraction of these negatives are therefore false, which pushes the measured
discrimination DOWN. The within-family figure is a lower bound, not a point estimate.
"""

from __future__ import annotations

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import numpy as np

from pipeline import config, seqtools
from pipeline.db import connect
from pipeline.embed import esm
from pipeline.negatives import sanity
from pipeline.train import dataset, head

MARKS = ",".join("?" * len(config.MEASURED_TIERS))


def fetch():
    with connect() as c:
        pos = list(c.execute(
            f"SELECT enzyme_id, sequence FROM characterised_enzymes WHERE is_positive=1 "
            f"AND excluded_from_training=0 AND sequence IS NOT NULL "
            f"AND family!='mhetase_like' AND source_ref IN ({MARKS})",
            config.MEASURED_TIERS))
        out_of_family = list(c.execute(
            "SELECT enzyme_id, sequence FROM characterised_enzymes "
            "WHERE is_negative=1 AND excluded_from_training=0 AND sequence IS NOT NULL"))
        near_miss = list(c.execute(
            "SELECT enzyme_id, sequence FROM characterised_enzymes "
            "WHERE is_near_miss=1 AND source_ref!='PAZy-nonPET' "
            "AND excluded_from_training=0 AND sequence IS NOT NULL"))
        within = list(c.execute(
            "SELECT enzyme_id, sequence FROM characterised_enzymes "
            "WHERE source_ref='PAZy-nonPET' AND sequence IS NOT NULL"))
        # The strongest negatives available: expressed, assayed under a published protocol,
        # and no product released. Kept as their OWN contrast rather than merged into the
        # PAZy set, because "measured inactive" and "not reported active" are different
        # claims and pooling them would let the weaker one ride on the stronger.
        measured_inactive = list(c.execute(
            "SELECT enzyme_id, sequence FROM characterised_enzymes "
            "WHERE within_family_basis='measured-inactive' AND excluded_from_training=0 "
            "AND sequence IS NOT NULL"))
        # Same negatives, sliced by which family definition admitted them. The conclusion
        # should not depend on that choice; if it does, that is worth knowing.
        by_basis = {}
        for b in ("cluster", "profile", "both"):
            q = ("SELECT enzyme_id, sequence FROM characterised_enzymes "
                 "WHERE source_ref='PAZy-nonPET' AND sequence IS NOT NULL AND "
                 + ("within_family_basis IN ('cluster','both')" if b == "cluster" else
                    "within_family_basis IN ('profile','both')" if b == "profile" else
                    "within_family_basis='both'"))
            by_basis[b] = [tuple(r) for r in c.execute(q)]
    return ([tuple(r) for r in pos], [tuple(r) for r in out_of_family],
            [tuple(r) for r in near_miss], [tuple(r) for r in within], by_basis,
            [tuple(r) for r in measured_inactive])


def evaluate(name, pos, neg, X_by_id, restrict_to_mixed=False):
    """Cluster-grouped evaluation of one positive/negative regime."""
    recs = list(pos) + list(neg)
    fa = seqtools.write_fasta(recs, config.INTERIM_DIR / f"wf_{name}.fasta")
    clu = seqtools.cluster(fa, min_seq_id=0.3)
    groups = np.array([clu.get(i, i) for i, _ in recs])

    if restrict_to_mixed:
        neg_clusters = {clu.get(i, i) for i, _ in neg}
        keep = [k for k, (i, _) in enumerate(recs) if clu.get(i, i) in neg_clusters]
        recs = [recs[k] for k in keep]
        y_all = np.array([1] * len(pos) + [0] * len(neg))[keep]
        groups = groups[keep]
    else:
        y_all = np.array([1] * len(pos) + [0] * len(neg))

    have = [k for k, (i, _) in enumerate(recs) if i in X_by_id]
    recs = [recs[k] for k in have]
    y = y_all[have]
    groups = groups[have]
    if len(set(y)) < 2:
        return {"regime": name, "note": "one class only after filtering"}

    Xt = np.vstack([X_by_id[i] for i, _ in recs])
    res, _ = head.train_eval(Xt, y, groups, name)

    base = sanity.run([s for (i, s), yy in zip(recs, y) if yy == 1],
                      [s for (i, s), yy in zip(recs, y) if yy == 0],
                      groups=list(groups))
    return {
        "regime": name,
        "n_pos": int((y == 1).sum()), "n_neg": int((y == 0).sum()),
        "n_clusters": len(set(groups)),
        "n_pos_clusters": len(set(groups[y == 1])),
        "n_neg_clusters": len(set(groups[y == 0])),
        "shared_clusters": len(set(groups[y == 1]) & set(groups[y == 0])),
        "folds": res.folds,
        "auc": round(res.auc, 4), "auc_std": round(res.auc_std, 4),
        "average_precision": round(res.average_precision, 4),
        "brier": round(res.brier, 4),
        "composition_auc": round(float(base["auc_mean"]), 4),
        "composition_null": round(float(base["null_auc_mean"]), 4),
    }


if __name__ == "__main__":
    dataset.apply_filters()
    pos, out_of_family, near_miss, within, by_basis, measured_inactive = fetch()
    print(f"measured positives        {len(pos)}")
    print(f"out-of-family negatives   {len(out_of_family)}")
    print(f"near-miss negatives       {len(near_miss)}")
    print(f"within-family, inferred   {len(within)}   (PAZy 'not reported active on PET')")
    print(f"within-family, MEASURED   {len(measured_inactive)}   "
          f"(expressed, assayed, no product released)\n", flush=True)

    ids, X = esm.load(config.PROCESSED_DIR / "embeddings.npz")
    X_by_id = {i: X[k] for k, i in enumerate(ids)}

    regimes = [
        ("out-of-family", out_of_family, False),
        ("near-miss", near_miss, False),
        ("within-family", within, False),
        ("within-family[MEASURED inactive]", measured_inactive, False),
        ("within-family[MEASURED, mixed clusters only]", measured_inactive, True),
        ("within-family-mixed-clusters-only", within, True),
        ("within-family[cluster-defined]", by_basis["cluster"], True),
        ("within-family[profile-defined]", by_basis["profile"], True),
        ("within-family[both-tests]", by_basis["both"], True),
    ]
    results = []
    for name, neg, restrict in regimes:
        print(f"--- {name} ---", flush=True)
        try:
            r = evaluate(name, pos, neg, X_by_id, restrict_to_mixed=restrict)
        except ValueError as exc:
            # Restricting to mixed clusters can leave a fold with one class in it. That is
            # a statement about how few clusters hold both, which is the finding this whole
            # script exists to measure -- so it is reported, not raised.
            print(json.dumps({"regime": name, "not_evaluable": str(exc)}, indent=2), flush=True)
            continue
        results.append(r)
        print(json.dumps(r, indent=2), flush=True)

    print("\n" + "=" * 104)
    print(f"{'regime':<36} {'pos':>5} {'neg':>5} {'shared':>7} "
          f"{'head AUC':>16} {'AP':>7} {'composition':>12}")
    print("=" * 104)
    for r in results:
        if "auc" not in r:
            continue
        print(f"{r['regime']:<36} {r['n_pos']:>5} {r['n_neg']:>5} {r['shared_clusters']:>7} "
              f"{r['auc']:>8.3f}+/-{r['auc_std']:<5.3f} {r['average_precision']:>7.3f} "
              f"{r['composition_auc']:>12.3f}")

    out = config.ROOT_DIR / "release" / "within_family_eval.json"
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {out}")
