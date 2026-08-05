import pathlib
import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import numpy as np
from scipy.stats import spearmanr
from pipeline import config, seqtools
from pipeline.db import connect
from pipeline.embed import esm
from pipeline.train import head

ids, X = esm.load(config.PROCESSED_DIR/"embeddings.npz")
idx = {i:k for k,i in enumerate(ids)}

with connect() as c:
    rows = list(c.execute(
        "SELECT enzyme_id, sequence, is_positive, is_negative, is_near_miss, source_ref "
        "FROM characterised_enzymes WHERE excluded_from_training=0 AND sequence IS NOT NULL "
        "AND family != 'mhetase_like'"))
rows = [r for r in rows if r["enzyme_id"] in idx]

eids = [r["enzyme_id"] for r in rows]
Xt   = np.vstack([X[idx[e]] for e in eids])
pos  = np.array([1 if r["is_positive"] else 0 for r in rows])
neg  = np.array([1 if (r["is_negative"] or r["is_near_miss"]) else 0 for r in rows])
evid = np.array([1 if (r["source_ref"] in head.EVIDENCED_TIERS and r["is_positive"]) else 0 for r in rows])

# Cluster at 30% identity: splits are by cluster, never by sequence (spec section 8).
fa = seqtools.write_fasta([(r["enzyme_id"], r["sequence"]) for r in rows], config.INTERIM_DIR/"train.fasta")
clu = seqtools.cluster(fa, min_seq_id=0.3)
groups = np.array([clu.get(e, e) for e in eids])

print(f"{len(rows)} sequences | positives {pos.sum()} (evidenced {evid.sum()}) | negatives {neg.sum()}")
print(f"clusters: {len(set(groups))} total, {len(set(groups[pos==1]))} carrying positives, "
      f"{len(set(groups[evid==1]))} carrying EVIDENCED positives\n")

results = {}

# --- naive: everything annotated positive is positive ---
mask = (pos == 1) | (neg == 1)
r, oof_naive = head.train_eval(Xt[mask], pos[mask], groups[mask], "naive")
results["naive"] = r

# --- evidence-only: auto-annotated dropped from training entirely ---
mask_e = (evid == 1) | (neg == 1)
r2, oof_ev = head.train_eval(Xt[mask_e], evid[mask_e], groups[mask_e], "evidence")
results["evidence"] = r2

# --- PU: auto-annotated treated as UNLABELLED, not positive ---
r3, oof_pu = head.train_eval_pu(Xt, evid, neg, groups)
results["pu"] = r3

print(f"{'scheme':<10} {'pos':>5} {'unlab':>6} {'neg':>5} {'clus':>5} {'folds':>6} {'AUC':>14} {'AP':>7} {'Brier':>7}  c")
for k,r in results.items():
    c = f"{r.c_estimate:.3f}" if r.c_estimate else "-"
    print(f"{k:<10} {r.n_pos:>5} {r.n_unlabelled:>6} {r.n_neg:>5} {r.n_pos_clusters:>5} "
          f"{r.folds:>6} {r.auc:>7.3f}+/-{r.auc_std:<5.3f} {r.average_precision:>7.3f} {r.brier:>7.3f}  {c}")

# --- do the schemes agree on how they rank the candidates? ---
print("\n=== candidate ranking agreement ===")
with connect() as c:
    cands = [r[0] for r in c.execute("SELECT candidate_id FROM candidates")]
cands = [c_ for c_ in cands if c_ in idx]
Xc = np.vstack([X[idx[c_]] for c_ in cands])

models = {}
models["naive"]    = head.fit_full(Xt[mask], pos[mask])
models["evidence"] = head.fit_full(Xt[mask_e], evid[mask_e])
scores = {k: m.predict_proba(Xc)[:,1] for k,m in models.items()}

rho, p = spearmanr(scores["naive"], scores["evidence"])
print(f"  naive vs evidence: Spearman rho = {rho:.3f} (p={p:.2g})")
for k,v in scores.items():
    print(f"  {k:<9} mean {v.mean():.3f}  median {np.median(v):.3f}  >0.5: {(v>0.5).sum()}/{len(v)}")

# top-10 overlap
top = {k: set(np.array(cands)[np.argsort(-v)[:10]]) for k,v in scores.items()}
print(f"  top-10 overlap: {len(top['naive'] & top['evidence'])}/10")

# --- baselines the head must beat ---
print("\n=== baselines ===")
with connect() as c:
    ev = {r[0]: r[1] for r in c.execute("SELECT candidate_id, recall_evalue FROM candidates WHERE recall_evalue IS NOT NULL")}
common = [c_ for c_ in cands if c_ in ev]
if common:
    e = np.array([-np.log10(max(ev[c_],1e-300)) for c_ in common])
    for k,v in scores.items():
        s = np.array([v[cands.index(c_)] for c_ in common])
        rr,_ = spearmanr(s, e)
        print(f"  {k:<9} vs E-value rank: Spearman rho = {rr:.3f}")
print("  composition baseline (cluster-grouped, from the gate): 0.778")
np.savez_compressed(config.PROCESSED_DIR/"candidate_scores.npz",
                    ids=np.array(cands,dtype=object),
                    naive=scores["naive"], evidence=scores["evidence"])
