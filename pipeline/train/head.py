"""The PET activity head, trained three ways so the labelling assumption is testable.

Spec section 5.3 warns that "not annotated as a PETase" is not "tested and inactive". The
converse bites just as hard here: **"annotated by similarity" is not "tested and active"**.
449 of the 529 positives carry EC 3.1.1.101 from `ECO:0000256` (automatic annotation), so
calling them positives asserts something no one measured.

If those are trained as confident positives, the head can score well by relearning "is
this in the polyesterase family", which is exactly what the annotation pipeline already
encoded. It would then look good on every internal metric while knowing nothing new.

So the head is trained under three labelling schemes and compared:

  naive     all 500 positives are positive. What the annotation says.
  evidence  only the 16 experimentally evidenced (ECO:0000269) are positive; the
            auto-annotated are dropped from training entirely. Small but honest.
  pu        Elkan-Noto: the 16 are labelled, the auto-annotated are UNLABELLED rather than
            positive, and P(y=1|x) = g(x)/c with c estimated on held-out labelled positives.

Agreement between them means the labels were fine. Divergence is the most informative
result available at this stage, and says whether more curation is worth the cost.

Everything is a shallow head on frozen embeddings. No fine-tuning: with order 10^1 real
positives, end-to-end training memorises (spec section 5.5).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .. import config

# Tiers whose label rests on an experiment rather than on similarity.
#
# THE THIRD COPY. This set was written before PAZy was integrated and never updated, so it
# omitted `PAZy-measured` — the tier carrying 320 of the 341 measured positives. The effect
# was not a wrong number but a silently useless one: run_train.py's evidence-only and PU
# schemes saw 13 evidenced positives in a single 30% cluster and recorded "NOT EVALUABLE",
# while the measured-only head scored 0.976 from a separate script with its own correct
# list. Two definitions of "evidenced", two different answers, no error anywhere.
#
# Now imported from pipeline.config, which derives it from VARIANTS. See the note there.
EVIDENCED_TIERS = set(config.MEASURED_TIERS)


@dataclass
class Result:
    scheme: str
    n_pos: int
    n_neg: int
    n_unlabelled: int
    n_pos_clusters: int
    auc: float
    auc_std: float
    average_precision: float
    brier: float
    c_estimate: Optional[float] = None
    folds: int = 0
    notes: str = ""


def _make_pipeline():
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler
    return make_pipeline(
        StandardScaler(),
        LogisticRegression(max_iter=3000, class_weight="balanced", C=1.0),
    )


def _grouped_splits(y: np.ndarray, groups: np.ndarray, n_splits: int, seed: int):
    from sklearn.model_selection import StratifiedGroupKFold
    n_pos_groups = len(set(groups[y == 1]))
    n_splits = max(2, min(n_splits, n_pos_groups))
    cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
    return list(cv.split(np.zeros(len(y)), y, groups=groups)), n_splits


def train_eval(X: np.ndarray, y: np.ndarray, groups: np.ndarray, scheme: str,
               n_splits: int = 4, seed: int = 0) -> Tuple[Result, np.ndarray]:
    """Cluster-grouped CV for the naive and evidence schemes.

    Returns the result and out-of-fold probabilities (NaN where a row was never in a
    test fold, which happens when its cluster carries no positives).
    """
    from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

    splits, n_splits = _grouped_splits(y, groups, n_splits, seed)
    oof = np.full(len(y), np.nan)
    aucs: List[float] = []

    for train, test in splits:
        if len(set(y[train])) < 2 or len(set(y[test])) < 2:
            continue
        clf = _make_pipeline()
        clf.fit(X[train], y[train])
        p = clf.predict_proba(X[test])[:, 1]
        oof[test] = p
        aucs.append(roc_auc_score(y[test], p))

    scored = ~np.isnan(oof)
    return Result(
        scheme=scheme,
        n_pos=int((y == 1).sum()), n_neg=int((y == 0).sum()), n_unlabelled=0,
        n_pos_clusters=len(set(groups[y == 1])),
        auc=float(np.mean(aucs)) if aucs else float("nan"),
        auc_std=float(np.std(aucs)) if aucs else float("nan"),
        average_precision=float(average_precision_score(y[scored], oof[scored])) if scored.any() else float("nan"),
        brier=float(brier_score_loss(y[scored], oof[scored])) if scored.any() else float("nan"),
        folds=len(aucs),
    ), oof


def train_eval_pu(X: np.ndarray, s: np.ndarray, y_true_neg: np.ndarray,
                  groups: np.ndarray, n_splits: int = 4, seed: int = 0
                  ) -> Tuple[Result, np.ndarray]:
    """Elkan-Noto PU.

    `s` is 1 for LABELLED positives and 0 for everything else (unlabelled plus known
    negatives). `y_true_neg` marks rows known to be genuine negatives, used only for
    honest evaluation, never for fitting.

    The non-traditional classifier g(x) predicts s. Under the selected-completely-at-random
    assumption, g(x) = c * P(y=1|x) where c = P(s=1|y=1), so c is estimated as the mean of
    g over held-out LABELLED positives and the calibrated estimate is g(x)/c.

    With 16 labelled positives that estimate is genuinely noisy, which is reported rather
    than hidden: c_estimate carries the number so it can be judged.
    """
    from sklearn.metrics import average_precision_score, brier_score_loss, roc_auc_score

    splits, n_splits = _grouped_splits(s, groups, n_splits, seed)
    oof = np.full(len(s), np.nan)
    cs: List[float] = []
    aucs: List[float] = []

    for train, test in splits:
        if len(set(s[train])) < 2:
            continue
        clf = _make_pipeline()
        clf.fit(X[train], s[train])

        held_pos = test[s[test] == 1]
        if len(held_pos) == 0:
            continue
        c = float(np.mean(clf.predict_proba(X[held_pos])[:, 1]))
        if c <= 1e-6:
            continue
        cs.append(c)

        p = np.clip(clf.predict_proba(X[test])[:, 1] / c, 0.0, 1.0)
        oof[test] = p

        # Evaluate on rows whose truth is actually known: labelled positives and known
        # negatives. The unlabelled have no ground truth, so scoring against them would
        # be measuring agreement with an assumption.
        known = test[(s[test] == 1) | (y_true_neg[test] == 1)]
        if len(known) and len(set(s[known])) > 1:
            aucs.append(roc_auc_score(s[known], np.clip(
                clf.predict_proba(X[known])[:, 1] / c, 0, 1)))

    scored = ~np.isnan(oof) & ((s == 1) | (y_true_neg == 1))
    return Result(
        scheme="pu",
        n_pos=int((s == 1).sum()), n_neg=int(y_true_neg.sum()),
        n_unlabelled=int(((s == 0) & (y_true_neg == 0)).sum()),
        n_pos_clusters=len(set(groups[s == 1])),
        auc=float(np.mean(aucs)) if aucs else float("nan"),
        auc_std=float(np.std(aucs)) if aucs else float("nan"),
        average_precision=float(average_precision_score(s[scored], oof[scored])) if scored.any() else float("nan"),
        brier=float(brier_score_loss(s[scored], oof[scored])) if scored.any() else float("nan"),
        c_estimate=float(np.mean(cs)) if cs else None,
        folds=len(aucs),
        notes=f"c estimated on {len(cs)} folds from {int((s==1).sum())} labelled positives",
    ), oof


def fit_full(X: np.ndarray, y: np.ndarray, calibrate: bool = True, seed: int = 0):
    """Fit on everything, with calibration, for scoring the candidates."""
    from sklearn.calibration import CalibratedClassifierCV

    base = _make_pipeline()
    if not calibrate:
        base.fit(X, y)
        return base
    n_pos = int((y == 1).sum())
    # Calibration needs enough positives per fold to mean anything; below that, isotonic
    # overfits badly and sigmoid is the safer choice.
    method = "isotonic" if n_pos >= 50 else "sigmoid"
    cv = max(2, min(3, n_pos))
    clf = CalibratedClassifierCV(base, method=method, cv=cv)
    clf.fit(X, y)
    return clf
