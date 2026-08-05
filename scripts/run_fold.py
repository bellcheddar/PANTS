import pathlib
import sys, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from pipeline import config
from pipeline.db import connect
from pipeline.db.manifest import stage_manifest
from pipeline.structure import fold

with connect() as c:
    cands = [(r[0], r[1]) for r in c.execute(
        "SELECT candidate_id, sequence FROM candidates ORDER BY recall_bitscore DESC")]
print(f"{len(cands)} candidates queued (cap {fold.MAX_FOLD_LENGTH} aa)", flush=True)
t0 = time.time()
with stage_manifest("structure", label="esmfold-v1-capped450",
                    model_version=config.ESMFOLD_MODEL,
                    params={"max_fold_length": fold.MAX_FOLD_LENGTH}) as m:
    res = fold.run(cands)
    m.counts(n_input=len(cands), n_output=len(res), n_discarded=len(cands)-len(res))
print(f"\nTOTAL {time.time()-t0:.0f}s for {len(res)} new structures", flush=True)
with connect() as c:
    print("structures:", c.execute("SELECT COUNT(*) FROM structures").fetchone()[0], flush=True)
    print("geometry:  ", c.execute("SELECT COUNT(*) FROM geometry").fetchone()[0], flush=True)
