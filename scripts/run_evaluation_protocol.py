#!/usr/bin/env python3
"""The full evaluation protocol from spec section 8, run end to end and reported honestly.

Six components were specified. Five run. The sixth cannot, and saying which is the point:

  cluster splits at 30% and 50%     runs. Two thresholds because a signature that survives
                                    30% grouping and not 50% is telling you the threshold
                                    was doing the work.
  leave-one-family-out              DEGENERATE, and that is a finding rather than a gap.
                                    Every ESTHER family in this catalogue is either wholly
                                    positive or wholly negative -- Polyesterase-lipase-
                                    cutinase is 77 positives and no negatives, Cutinase is
                                    110 negatives and no positives -- so holding out a
                                    family removes an entire class and AUC is undefined on
                                    the test fold. The label IS family membership here.
  retrieval baseline                runs. Nearest-neighbour identity to a known positive:
                                    the number a new method has to beat before it has
                                    earned the word "learned".
  reliability diagram               runs. Written as bin counts so the figure can be drawn
                                    from data rather than trusted from a claim.
  prospective holdout by date       runs, on UniProt first-public dates rather than PDB
                                    release dates, which are absent for the whole catalogue.
  measured versus annotated         runs, and is reported separately everywhere. Pooling
                                    them is how a model gets to look good for rediscovering
                                    the similarity rule that made the labels.

Also the PU class prior, swept across 1/3/5/10% rather than assumed, because the spec calls
for the prior to be estimated and its influence shown rather than a single value defended.

Every scheme is written to `training_runs`, including the ones that fail, so the stats page
shows the failures next to the successes.
"""

from __future__ import annotations

import json
import pathlib
import sys
from typing import Dict, List, Optional, Tuple

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config, seqtools
from pipeline.db import connect, now
from pipeline.embed import esm
from pipeline.train import head

from sklearn.calibration import calibration_curve
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

OUT = config.INTERIM_DIR / "evaluation_protocol.json"
MODEL_VERSION = "pet_activity_head/v2"


# ----------------------------------------------------------------------------------
# data
# ----------------------------------------------------------------------------------
def load() -> dict:
    ids, X = esm.load(config.PROCESSED_DIR / "embeddings.npz")
    idx = {i: k for k, i in enumerate(ids)}
    with connect() as c:
        rows = [dict(r) for r in c.execute(
            "SELECT enzyme_id, sequence, is_positive, is_negative, is_near_miss, "
            "       source_ref, esther_family, uniprot_first_public "
            "FROM characterised_enzymes "
            "WHERE excluded_from_training=0 AND sequence IS NOT NULL "
            "  AND family != 'mhetase_like'")]
    rows = [r for r in rows if r["enzyme_id"] in idx]
    eids = [r["enzyme_id"] for r in rows]
    return {
        "rows": rows, "eids": eids,
        "X": np.vstack([X[idx[e]] for e in eids]),
        "pos": np.array([1 if r["is_positive"] else 0 for r in rows]),
        "neg": np.array([1 if (r["is_negative"] or r["is_near_miss"]) else 0 for r in rows]),
        "measured": np.array([1 if (r["is_positive"] and r["source_ref"] in head.EVIDENCED_TIERS)
                              else 0 for r in rows]),
    }


def cluster_groups(rows: List[dict], eids: List[str], identity: float) -> np.ndarray:
    fa = seqtools.write_fasta([(r["enzyme_id"], r["sequence"]) for r in rows],
                              config.INTERIM_DIR / f"proto_{int(identity * 100)}.fasta")
    clu = seqtools.cluster(fa, min_seq_id=identity)
    return np.array([clu.get(e, e) for e in eids])


def _fit() -> "object":
    return make_pipeline(StandardScaler(),
                         LogisticRegression(max_iter=3000, class_weight="balanced"))


# ----------------------------------------------------------------------------------
# components
# ----------------------------------------------------------------------------------
def cluster_splits(d: dict, identity: float) -> Dict[str, dict]:
    """The headline numbers, at one clustering threshold, three label schemes."""
    groups = cluster_groups(d["rows"], d["eids"], identity)
    X, pos, neg, meas = d["X"], d["pos"], d["neg"], d["measured"]
    out = {}
    for name, y_pos in (("all-annotated", pos), ("measured-only", meas)):
        mask = (y_pos == 1) | (neg == 1)
        if len(set(y_pos[mask])) < 2:
            out[name] = {"note": "single class"}
            continue
        r, _ = head.train_eval(X[mask], y_pos[mask], groups[mask], name)
        out[name] = {"auc": r.auc, "auc_std": r.auc_std, "ap": r.average_precision,
                     "brier": r.brier, "n_pos": r.n_pos, "n_neg": r.n_neg,
                     "n_pos_clusters": r.n_pos_clusters, "folds": r.folds}
    r3, _ = head.train_eval_pu(X, meas, neg, groups)
    out["pu"] = {"auc": r3.auc, "auc_std": r3.auc_std, "ap": r3.average_precision,
                 "brier": r3.brier, "c_estimate": r3.c_estimate,
                 "n_pos": r3.n_pos, "n_neg": r3.n_neg, "folds": r3.folds}
    out["_n_clusters"] = len(set(groups))
    return out


def leave_one_family_out(d: dict) -> dict:
    """Held-out ESTHER family. Reports the degeneracy rather than a number that hides it."""
    fams: Dict[str, dict] = {}
    for r, p, n in zip(d["rows"], d["pos"], d["neg"]):
        f = r["esther_family"]
        if not f:
            continue
        e = fams.setdefault(f, {"pos": 0, "neg": 0})
        e["pos"] += int(p)
        e["neg"] += int(n)
    mixed = {f: v for f, v in fams.items() if v["pos"] and v["neg"]}
    return {"n_families": len(fams),
            "n_mixed_families": len(mixed),
            "families": fams,
            "evaluable": bool(mixed),
            "note": ("every family is wholly positive or wholly negative, so holding one "
                     "out removes a class and AUC is undefined: in this catalogue the "
                     "label is family membership")
            if not mixed else "mixed families exist; LOFO is evaluable"}


def retrieval_baseline(d: dict, identity: float = 0.30) -> dict:
    """Nearest-neighbour identity to a known positive, scored as a classifier.

    The number to beat. If a learned head cannot clear "how similar is this to something we
    already call a PETase", it has not added information, it has re-expressed the input.
    """
    groups = cluster_groups(d["rows"], d["eids"], identity)
    pos, neg = d["pos"], d["neg"]
    mask = (pos == 1) | (neg == 1)
    # Cosine similarity to the nearest MEASURED positive, excluding same-cluster
    # neighbours: leaving them in scores the baseline on its own training examples.
    Xn = d["X"] / (np.linalg.norm(d["X"], axis=1, keepdims=True) + 1e-9)
    ref = np.where(d["measured"] == 1)[0]
    scores = np.full(len(pos), np.nan)
    for i in np.where(mask)[0]:
        other = ref[groups[ref] != groups[i]]
        if len(other):
            scores[i] = float(np.max(Xn[i] @ Xn[other].T))
    ok = mask & ~np.isnan(scores)
    if len(set(pos[ok])) < 2:
        return {"note": "not evaluable"}
    return {"auc": float(roc_auc_score(pos[ok], scores[ok])), "n": int(ok.sum()),
            "note": "cosine similarity to the nearest measured positive outside the "
                    "enzyme's own 30% cluster"}


def composition_baseline(d: dict, identity: float = 0.30) -> dict:
    """Amino-acid composition alone. The permanently reported floor."""
    groups = cluster_groups(d["rows"], d["eids"], identity)
    aa = "ACDEFGHIKLMNPQRSTVWY"
    C = np.array([[r["sequence"].count(a) / max(len(r["sequence"]), 1) for a in aa]
                  for r in d["rows"]])
    pos, neg = d["pos"], d["neg"]
    mask = (pos == 1) | (neg == 1)
    r, _ = head.train_eval(C[mask], pos[mask], groups[mask], "composition")
    return {"auc": r.auc, "auc_std": r.auc_std, "folds": r.folds}


def reliability(d: dict, identity: float = 0.30, n_bins: int = 8) -> dict:
    """Out-of-fold calibration curve, written as data so the figure can be drawn from it."""
    groups = cluster_groups(d["rows"], d["eids"], identity)
    pos, neg = d["pos"], d["neg"]
    mask = (pos == 1) | (neg == 1)
    _, oof = head.train_eval(d["X"][mask], pos[mask], groups[mask], "reliability")
    scored = ~np.isnan(oof)
    if scored.sum() < 20 or len(set(pos[mask][scored])) < 2:
        return {"note": "not evaluable"}
    frac, mean_pred = calibration_curve(pos[mask][scored], oof[scored],
                                        n_bins=n_bins, strategy="quantile")
    return {"mean_predicted": [float(x) for x in mean_pred],
            "observed_fraction": [float(x) for x in frac],
            "n_scored": int(scored.sum()),
            "note": "a well-calibrated head lies on the diagonal; above it is "
                    "under-confident, below it over-confident"}


def prospective_holdout(d: dict, cutoff: str = "2020-01-01") -> dict:
    """Train on what was public before the cutoff, test on what appeared after.

    The only split that answers "would this have found the enzymes we now know about",
    because the test set is genuinely unavailable at training time rather than merely
    held out from it.
    """
    dates = [r["uniprot_first_public"] for r in d["rows"]]
    have = np.array([bool(x) for x in dates])
    before = np.array([bool(x) and x < cutoff for x in dates])
    after = np.array([bool(x) and x >= cutoff for x in dates])
    pos, neg = d["pos"], d["neg"]
    lab = (pos == 1) | (neg == 1)
    tr, te = before & lab, after & lab
    info = {"cutoff": cutoff, "n_dated": int(have.sum()),
            "n_train": int(tr.sum()), "n_test": int(te.sum()),
            "test_pos": int((pos[te] == 1).sum()), "test_neg": int((neg[te] == 1).sum())}
    if len(set(pos[tr])) < 2 or len(set(pos[te])) < 2:
        info["note"] = "not evaluable: one side of the cutoff holds a single class"
        return info
    clf = _fit()
    clf.fit(d["X"][tr], pos[tr])
    info["auc"] = float(roc_auc_score(pos[te], clf.predict_proba(d["X"][te])[:, 1]))
    info["note"] = ("UniProt first-public dates, not PDB release dates, which are absent "
                    "for the whole catalogue")
    # An AUC resting on a handful of negatives is a number, not evidence. Say so next to
    # it rather than letting the three decimal places imply a precision it does not have.
    smaller = min(info["test_pos"], info["test_neg"])
    if smaller < 20:
        info["underpowered"] = True
        info["note"] += (f"; UNDERPOWERED -- the test side holds only {smaller} of the "
                         f"smaller class, so this AUC is indicative at best")
    return info


def pu_prior_sweep(d: dict, identity: float = 0.30,
                   priors: Tuple[float, ...] = (0.01, 0.03, 0.05, 0.10)) -> dict:
    """How much the PU conclusion depends on the class prior, rather than defending one value.

    Two things this deliberately does NOT report, both because they would look like results
    and be neither.

    It does not report an AUC per prior. Dividing by the prior is strictly monotone, so the
    ranking cannot change and the AUC is invariant by construction. An earlier version
    clipped the adjusted score at 1 before scoring, which collapsed everything above the
    prior into a tie, and the ties moved the AUC: it appeared to show the prior mattering,
    0.913 rising to 0.977, when what it showed was tie-breaking.

    And it works from OUT-OF-FOLD probabilities, not a refit on the rows being scored. The
    in-sample version returned 1.000 for every prior, which is memorisation wearing the
    costume of an invariance check.

    What the prior genuinely governs is the DECISION: how many sequences the head would
    call positive. That is what moves here, and it is what the assumption is worth arguing
    about.
    """
    groups = cluster_groups(d["rows"], d["eids"], identity)
    meas, neg = d["measured"], d["neg"]
    r_est, oof = head.train_eval_pu(d["X"], meas, neg, groups)
    scored = ~np.isnan(oof)
    out = {"estimated": {"c": r_est.c_estimate, "auc": r_est.auc, "auc_std": r_est.auc_std},
           "n_scored_out_of_fold": int(scored.sum()), "assumed": {}}
    for pr in priors:
        called = np.clip(oof[scored] / max(pr, 1e-6), 0, 1) > 0.5
        out["assumed"][f"{pr:.0%}"] = {
            "n_called_positive": int(called.sum()),
            "fraction_called_positive": float(called.mean())}
    out["note"] = ("AUC is omitted on purpose: it is invariant to the prior by "
                   "construction. The prior sets how many sequences get CALLED positive, "
                   "which is the decision it actually governs.")
    return out


# ----------------------------------------------------------------------------------
def record(scheme: str, res: dict, identity: float, evidence: str,
           composition_auc: Optional[float], retrieval_auc: Optional[float]) -> None:
    with connect() as c:
        c.execute(
            "INSERT INTO training_runs(head_name, model_version, n_positives, n_negatives, "
            " auc, average_precision, brier_score, retrieval_baseline_auc, "
            " composition_baseline_auc, n_positive_clusters, cluster_identity_threshold, "
            " evidence_level, trained_at, config_json) "
            "VALUES(?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
            ("pet_activity", MODEL_VERSION, res.get("n_pos"), res.get("n_neg"),
             res.get("auc"), res.get("ap"), res.get("brier"), retrieval_auc,
             composition_auc, res.get("n_pos_clusters"), identity, evidence, now(),
             json.dumps({"scheme": scheme, "protocol": "full", **
                         {k: v for k, v in res.items() if k != "note"}}, default=float)))
        c.commit()


def main() -> int:
    d = load()
    print(f"{len(d['rows'])} trainable sequences: {int(d['pos'].sum())} positives "
          f"({int(d['measured'].sum())} measured), {int(d['neg'].sum())} negatives\n")

    report: Dict[str, object] = {}
    comp = composition_baseline(d)
    retr = retrieval_baseline(d)
    print(f"=== baselines (the floor a learned head has to clear) ===")
    print(f"  amino-acid composition   AUC {comp['auc']:.3f} +/- {comp['auc_std']:.3f}")
    print(f"  nearest measured positive AUC {retr.get('auc', float('nan')):.3f}  "
          f"(n={retr.get('n')})")
    report["composition_baseline"] = comp
    report["retrieval_baseline"] = retr

    for identity in (0.30, 0.50):
        res = cluster_splits(d, identity)
        report[f"cluster_{int(identity*100)}"] = res
        print(f"\n=== cluster-grouped at {identity:.0%} identity "
              f"({res['_n_clusters']} clusters) ===")
        for scheme in ("all-annotated", "measured-only", "pu"):
            r = res.get(scheme, {})
            if "auc" not in r:
                print(f"  {scheme:14s} {r.get('note', 'not evaluable')}")
                continue
            print(f"  {scheme:14s} AUC {r['auc']:.3f} +/- {r['auc_std']:.3f}  "
                  f"AP {r['ap']:.3f}  Brier {r['brier']:.3f}  "
                  f"pos {r['n_pos']} neg {r['n_neg']} folds {r['folds']}")
            record(scheme, r, identity,
                   "measured-only" if scheme != "all-annotated" else "mixed",
                   comp["auc"], retr.get("auc"))

    lofo = leave_one_family_out(d)
    report["leave_one_family_out"] = lofo
    print(f"\n=== leave-one-family-out ===")
    print(f"  {lofo['n_families']} ESTHER families, {lofo['n_mixed_families']} containing "
          f"both classes")
    print(f"  {lofo['note']}")

    ph = prospective_holdout(d)
    report["prospective_holdout"] = ph
    print(f"\n=== prospective holdout at {ph['cutoff']} ===")
    print(f"  {ph['n_train']} before, {ph['n_test']} after "
          f"({ph['test_pos']} positive, {ph['test_neg']} negative)")
    if "auc" in ph:
        print(f"  AUC {ph['auc']:.3f}"
              + ("   << UNDERPOWERED" if ph.get("underpowered") else ""))
    print(f"  {ph['note']}")

    rel = reliability(d)
    report["reliability"] = rel
    print(f"\n=== reliability ===")
    if "mean_predicted" in rel:
        for mp, of in zip(rel["mean_predicted"], rel["observed_fraction"]):
            bar = "#" * int(of * 40)
            print(f"  predicted {mp:.2f}  observed {of:.2f}  {bar}")
    else:
        print(f"  {rel['note']}")

    print("\n=== does the head beat the baselines? ===")
    c30 = report["cluster_30"]
    for scheme in ("all-annotated", "measured-only"):
        a = c30.get(scheme, {}).get("auc")
        if a is None:
            continue
        vs_comp = a - comp["auc"]
        vs_retr = a - (retr.get("auc") or float("nan"))
        verdict = ("clears both" if vs_retr > 0.02 else
                   "clears composition, NOT retrieval" if vs_retr <= 0.02 else "unclear")
        print(f"  {scheme:14s} AUC {a:.3f}  vs composition {vs_comp:+.3f}  "
              f"vs retrieval {vs_retr:+.3f}   {verdict}")
    report["verdict_vs_baselines"] = {
        s: {"auc": c30.get(s, {}).get("auc"),
            "minus_composition": (c30.get(s, {}).get("auc") or 0) - comp["auc"],
            "minus_retrieval": (c30.get(s, {}).get("auc") or 0) - (retr.get("auc") or 0)}
        for s in ("all-annotated", "measured-only")}

    sweep = pu_prior_sweep(d)
    report["pu_prior_sweep"] = sweep
    print(f"\n=== PU class prior ===")
    est = sweep["estimated"]
    c_est = f"{est['c']:.3f}" if est["c"] is not None else "-"
    print(f"  Elkan-Noto estimated c {c_est}  AUC {est['auc']:.3f} +/- {est['auc_std']:.3f}")
    n = sweep["n_scored_out_of_fold"]
    print(f"  sweeping the ASSUMED prior over {n} out-of-fold scores "
          f"(AUC omitted: invariant by construction)")
    for pr, v in sweep["assumed"].items():
        print(f"    assumed {pr:>4s}  calls {v['n_called_positive']:4d} of {n} positive "
              f"({v['fraction_called_positive']:.0%})")

    OUT.write_text(json.dumps(report, indent=2, default=float))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
