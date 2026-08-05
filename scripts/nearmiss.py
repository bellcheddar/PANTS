"""The harder test: measured PET-active enzymes vs NEAR MISSES.

The 0.976 result used hard negatives from other alpha/beta-hydrolase families, which is a
family-level question. This asks the one that matters: can the head separate a measured
PET degrader from a cutinase-family esterase that works on soluble esters but is not a
meaningful degrader of crystalline PET?
"""
import pathlib
import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import numpy as np
from pipeline import config, seqtools
from pipeline.db import connect
from pipeline.embed import esm
from pipeline.negatives import sanity
from pipeline.train import head, dataset

# Fifth copy of this list, now removed. See pipeline/config.py for why.
MEASURED = config.MEASURED_TIERS
marks = ','.join('?'*len(MEASURED))
dataset.apply_filters()

with connect() as c:
    pos = list(c.execute(
        f"SELECT enzyme_id, sequence FROM characterised_enzymes WHERE is_positive=1 "
        f"AND excluded_from_training=0 AND sequence IS NOT NULL AND family!='mhetase_like' "
        f"AND source_ref IN ({marks})", MEASURED))
    nm = list(c.execute(
        "SELECT enzyme_id, sequence FROM characterised_enzymes WHERE is_near_miss=1 "
        "AND excluded_from_training=0 AND sequence IS NOT NULL"))
    far = list(c.execute(
        "SELECT enzyme_id, sequence FROM characterised_enzymes WHERE is_negative=1 "
        "AND excluded_from_training=0 AND sequence IS NOT NULL"))

def run(name, positives, negatives):
    recs = [(r[0], r[1]) for r in positives] + [(r[0], r[1]) for r in negatives]
    fa = seqtools.write_fasta(recs, config.INTERIM_DIR / f"nm_{name}.fasta")
    clu = seqtools.cluster(fa, min_seq_id=0.3)
    groups = np.array([clu.get(i, i) for i, _ in recs])
    npc = len(set(groups[:len(positives)]))

    comp = sanity.run([r[1] for r in positives], [r[1] for r in negatives], groups=list(groups))

    ids, X = esm.load(config.PROCESSED_DIR / "embeddings.npz")
    have = {i: k for k, i in enumerate(ids)}
    keep = [n for n, (i, _) in enumerate(recs) if i in have]
    Xt = np.vstack([X[have[recs[n][0]]] for n in keep])
    y = np.array([1 if n < len(positives) else 0 for n in keep])
    res, _ = head.train_eval(Xt, y, groups[keep], name)
    print(f"\n=== {name} ===", flush=True)
    print(f"  {len(positives)} positives ({npc} clusters) vs {len(negatives)} negatives", flush=True)
    print(f"  ESM-2 head        AUC {res.auc:.3f} +/- {res.auc_std:.3f}   AP {res.average_precision:.3f}   Brier {res.brier:.3f}   folds {res.folds}", flush=True)
    print(f"  composition only  AUC {comp['auc_mean']:.3f} +/- {comp['auc_std']:.3f}   (null {comp['null_auc_mean']:.3f})", flush=True)
    print(f"  margin over composition: {res.auc - comp['auc_mean']:+.3f}", flush=True)

run("vs distant families (the 0.976 result)", pos, far)
run("vs NEAR MISSES (the hard question)", pos, nm)
run("vs both", pos, far + nm)
