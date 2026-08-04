"""Risk 1: the trivial-baseline check on the hard negative set.

PLAN_v1.md ranks this as the single point of failure for the whole project, and it is the
one check that must run BEFORE any ESM-2 work, because it is the only cheap way to find
out that the negative set is broken.

The logic: fit a deliberately stupid classifier (amino-acid composition and length, no
structure, no embedding, no fold information) on positives versus negatives. If that
already separates them well, the model does not need to learn anything about polyester
chemistry to score well, and every downstream number would look good while meaning
nothing. A high AUC here is a FAILURE, not a success.

With order 10^1 positives the estimate is noisy by construction, so this reports repeated
stratified cross-validation with a spread, and compares against a label-shuffled null.
A single point estimate at this sample size would be close to meaningless.
"""

from __future__ import annotations

from typing import Dict, List, Optional, Tuple

import numpy as np

AA = "ACDEFGHIKLMNPQRSTVWY"

# Above this, composition alone separates the classes well enough that the negative set
# cannot be trusted: the matching has failed and must be redone.
FAIL_THRESHOLD = 0.90
# Below this the sets are compositionally indistinguishable, which is the ideal.
PASS_THRESHOLD = 0.75


def featurise(sequences: List[str]) -> np.ndarray:
    """20 amino-acid fractions plus length. Nothing else, on purpose."""
    rows = []
    for s in sequences:
        n = max(1, len(s))
        rows.append([s.count(a) / n for a in AA] + [len(s)])
    return np.asarray(rows, dtype=float)


def run(pos_seqs: List[str], neg_seqs: List[str], n_repeats: int = 20,
        seed: int = 0, groups: Optional[List[str]] = None) -> Dict[str, object]:
    """Cross-validated trivial-baseline AUC.

    `groups` (one label per sequence, positives then negatives) makes the split
    cluster-aware. Without it, engineered variants of the same parent land in both train
    and test and the AUC is a leakage artefact rather than a measurement: the nine curated
    positives are 98%+ identical in places and form one cluster, so ungrouped CV on them
    reports ~1.0 no matter how good the negative set is.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.metrics import roc_auc_score
    from sklearn.model_selection import RepeatedStratifiedKFold, StratifiedGroupKFold
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    X = featurise(pos_seqs + neg_seqs)
    y = np.array([1] * len(pos_seqs) + [0] * len(neg_seqs))
    g = np.asarray(groups) if groups is not None else None

    n_pos_groups = len(set(g[y == 1])) if g is not None else len(pos_seqs)
    # n_splits is bounded by the number of INDEPENDENT positive units, not the raw count.
    n_splits = max(2, min(4, n_pos_groups // 2)) if n_pos_groups >= 4 else 2

    def cv_auc(labels: np.ndarray, rs: int) -> List[float]:
        aucs: List[float] = []
        if g is not None:
            cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=rs)
            splits = list(cv.split(X, labels, groups=g))
        else:
            cv = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats,
                                         random_state=rs)
            splits = list(cv.split(X, labels))
        for train, test in splits:
            if len(set(labels[test])) < 2:
                continue
            clf = make_pipeline(
                StandardScaler(),
                LogisticRegression(max_iter=2000, class_weight="balanced"),
            )
            clf.fit(X[train], labels[train])
            aucs.append(roc_auc_score(labels[test], clf.predict_proba(X[test])[:, 1]))
        return aucs

    real = cv_auc(y, seed)
    # Label-shuffled null: how high does this AUC go on pure noise at this sample size?
    rng = np.random.default_rng(seed)
    shuffled = y.copy()
    rng.shuffle(shuffled)
    null = cv_auc(shuffled, seed + 1)

    mean = float(np.mean(real))
    verdict = ("FAIL" if mean >= FAIL_THRESHOLD
               else "PASS" if mean <= PASS_THRESHOLD
               else "MARGINAL")

    return {
        "n_positives": len(pos_seqs),
        "n_negatives": len(neg_seqs),
        "n_splits": n_splits,
        "auc_mean": round(mean, 4),
        "auc_std": round(float(np.std(real)), 4),
        "auc_p10": round(float(np.percentile(real, 10)), 4),
        "auc_p90": round(float(np.percentile(real, 90)), 4),
        "null_auc_mean": round(float(np.mean(null)), 4),
        "null_auc_p90": round(float(np.percentile(null, 90)), 4),
        "verdict": verdict,
        "n_folds_scored": len(real),
    }


def top_discriminating_features(pos_seqs: List[str], neg_seqs: List[str],
                                k: int = 6) -> List[Tuple[str, float]]:
    """Which composition features carry the separation, as a standardised coefficient.

    Diagnostic, not a metric: if the separation is real, this says what to fix.
    """
    from sklearn.linear_model import LogisticRegression
    from sklearn.pipeline import make_pipeline
    from sklearn.preprocessing import StandardScaler

    X = featurise(pos_seqs + neg_seqs)
    y = np.array([1] * len(pos_seqs) + [0] * len(neg_seqs))
    clf = make_pipeline(StandardScaler(),
                        LogisticRegression(max_iter=2000, class_weight="balanced"))
    clf.fit(X, y)
    coefs = clf[-1].coef_[0]
    names = [f"{a} fraction" for a in AA] + ["length"]
    order = np.argsort(np.abs(coefs))[::-1][:k]
    return [(names[i], round(float(coefs[i]), 3)) for i in order]
