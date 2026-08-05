import pathlib
import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
import numpy as np
from pipeline import config
from pipeline.db import connect
from pipeline.db.manifest import stage_manifest
from pipeline.embed import esm
from pipeline.train import dataset

recs, labels, tiers = dataset.training_set()
with connect() as c:
    cands = [(r[0], r[1]) for r in c.execute("SELECT candidate_id, sequence FROM candidates")]
print(f"training {len(recs)}, candidates {len(cands)}", flush=True)

with stage_manifest("embed", label="esm2-t12-35M", model_version=config.ESM_MODEL) as m:
    all_recs = recs + cands
    ids, X, rep = esm.embed(all_recs, batch_size=8)
    print(f"\n{rep['n']} embedded, dim {rep['dim']}, {rep['seconds']}s "
          f"({rep['seq_per_sec']} seq/s), truncated {len(rep['truncated'])}", flush=True)

    out = config.PROCESSED_DIR / "embeddings.npz"
    esm.save(ids, X, out)

    idx = {i: k for k, i in enumerate(ids)}
    np.savez_compressed(config.PROCESSED_DIR / "train_labels.npz",
                        ids=np.array([r[0] for r in recs], dtype=object),
                        labels=np.array(labels), tiers=np.array(tiers, dtype=object))
    m.counts(n_input=len(all_recs), n_output=rep["n"], n_discarded=len(rep["truncated"]))
    print("saved", out, f"{out.stat().st_size/1e6:.1f} MB", flush=True)

# Sanity: do the embeddings put known relatives near each other?
from numpy.linalg import norm
Xn = X / np.clip(norm(X, axis=1, keepdims=True), 1e-9, None)
pos = {"IsPETase","FAST-PETase","LCC","LCC-ICCG","TfCut2","Cut190"}
have = [p for p in pos if p in idx]
print("\ncosine similarity among known enzymes:")
for a in have:
    sims = sorted(((float(Xn[idx[a]] @ Xn[idx[b]]), b) for b in have if b != a), reverse=True)
    print(f"  {a:<14} nearest: " + ", ".join(f"{b} {s:.3f}" for s, b in sims[:3]))
