import pathlib
import sys, pickle, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from pathlib import Path
from pipeline import config
from pipeline.db import connect
from pipeline.recall import run as recall_run, library

lib = pickle.load(open(config.INTERIM_DIR/"library.pkl","rb"))
with connect() as c:
    pos = [(r[0], r[1]) for r in c.execute(
        "SELECT enzyme_id, sequence FROM characterised_enzymes "
        "WHERE is_positive=1 AND sequence IS NOT NULL AND family!='mhetase_like'")]

ENV = {
 "ERZ782921":  "landfill",
 "ERZ21854033":"marine_plastisphere",
 "ERZ2185403": "marine_plastisphere",
 "ERZ794970":  "compost",
 "ERZ10545954":"compost",
}
files = sorted((config.RAW_DIR/"mgnify").glob("*.faa"))
print(f"{len(files)} FASTA files, {len(pos)} positives, {len(lib)} profiles", flush=True)

t0=time.time()
for f in files:
    env = next((v for k,v in ENV.items() if f.name.startswith(k)), "unknown")
    t=time.time()
    res = recall_run.run([f], lib, pos, environment=env, label="mgnify-v1")
    print(f"{f.name[:34]:<36} env={env:<20} scanned={res.n_scanned:>7,} "
          f"prefilter={res.n_prefilter:>5} profile={res.n_profile_matched:>5} "
          f"triad={res.n_triad_complete:>5} written={res.n_written:>5}  {time.time()-t:.0f}s", flush=True)
print(f"TOTAL {time.time()-t0:.0f}s")
with connect() as c:
    print("candidates in db:", c.execute("SELECT COUNT(*) FROM candidates").fetchone()[0])
    for r in c.execute("SELECT source_environment, COUNT(*) n FROM candidates GROUP BY 1 ORDER BY n DESC"):
        print(f"  {r[0]:<22} {r[1]}")
