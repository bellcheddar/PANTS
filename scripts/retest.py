import pathlib
import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import numpy as np
from pipeline import config, seqtools
from pipeline.db import connect
from pipeline.embed import esm
from pipeline.negatives import sanity
from pipeline.train import head, dataset

# Fourth copy of this list, now removed. See pipeline/config.py for why.
MEASURED = config.MEASURED_TIERS

dataset.apply_filters()
marks = ','.join('?'*len(MEASURED))
with connect() as c:
    pos = list(c.execute(
        f"SELECT enzyme_id, sequence FROM characterised_enzymes WHERE is_positive=1 "
        f"AND excluded_from_training=0 AND sequence IS NOT NULL AND family!='mhetase_like' "
        f"AND source_ref IN ({marks})", MEASURED))
    neg = list(c.execute(
        "SELECT enzyme_id, sequence FROM characterised_enzymes "
        "WHERE (is_negative=1 OR is_near_miss=1) AND excluded_from_training=0 "
        "AND sequence IS NOT NULL"))
print(f"measured positives {len(pos)} | negatives+near-misses {len(neg)}", flush=True)

recs = [(r[0], r[1]) for r in pos] + [(r[0], r[1]) for r in neg]
fa = seqtools.write_fasta(recs, config.INTERIM_DIR/"retest.fasta")
clu = seqtools.cluster(fa, min_seq_id=0.3)
groups = np.array([clu.get(i, i) for i, _ in recs])
print(f"positives span {len(set(groups[:len(pos)]))} clusters at 30%", flush=True)

# ---- risk 1: composition-only baseline, cluster-grouped ----
r = sanity.run([r[1] for r in pos], [r[1] for r in neg], groups=list(groups))
print("\n=== RISK 1: composition baseline, measured positives only ===", flush=True)
for k in ("auc_mean","auc_std","auc_p10","auc_p90","null_auc_mean","verdict"):
    print(f"  {k:<16}{r[k]}", flush=True)

# ---- the head, on evidence only, which was previously NOT EVALUABLE ----
ids, X = esm.load(config.PROCESSED_DIR/"embeddings.npz")
have = {i: k for k, i in enumerate(ids)}
missing = [i for i, _ in recs if i not in have]
if missing:
    print(f"\nembedding {len(missing)} new sequences...", flush=True)
    seqmap = dict(recs)
    nids, nX, rep = esm.embed([(m, seqmap[m]) for m in missing], batch_size=8, progress_every=0)
    ids = list(ids) + list(nids); X = np.vstack([X, nX])
    esm.save(ids, X, config.PROCESSED_DIR/"embeddings.npz")
    have = {i: k for k, i in enumerate(ids)}
    print(f"  done in {rep['seconds']}s", flush=True)

keep = [n for n, (i, _) in enumerate(recs) if i in have]
Xt = np.vstack([X[have[recs[n][0]]] for n in keep])
y  = np.array([1 if n < len(pos) else 0 for n in keep])
g  = groups[keep]
res, _ = head.train_eval(Xt, y, g, "evidence-only")
print("\n=== HEAD on measured positives (was 'not evaluable') ===", flush=True)
for k in ("n_pos","n_neg","n_pos_clusters","folds","auc","auc_std","average_precision","brier"):
    print(f"  {k:<20}{getattr(res,k)}", flush=True)
