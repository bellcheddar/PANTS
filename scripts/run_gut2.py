import pathlib
import sys, pickle, time, json
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from pathlib import Path
from pipeline import config
from pipeline.db import connect
from pipeline.recall import run as recall_run, mgnify

lib = pickle.load(open(config.INTERIM_DIR/"library2.pkl","rb"))
with connect() as c:
    pos = [(r[0], r[1]) for r in c.execute(
        "SELECT enzyme_id, sequence FROM characterised_enzymes "
        "WHERE is_positive=1 AND excluded_from_training=0 AND sequence IS NOT NULL "
        "AND family!='mhetase_like'")]

DONE = Path(config.INTERIM_DIR/"gut_done.txt")
done = set(DONE.read_text().split()) if DONE.exists() else set()
files = [p for p in sorted((config.RAW_DIR/"mgnify").glob("*.fa*"))
         if mgnify.environment_for(p.name) == "human_gut" and p.name not in done]
print(f"{len(files)} gut files to scan, {len(pos)} positives (HGMPs now included)", flush=True)
t0=time.time()
for i, f in enumerate(files, 1):
    t=time.time()
    res = recall_run.run([f], lib, pos, environment="human_gut", label="gut-v2")
    with DONE.open("a") as fh: fh.write(f.name+"\n")
    print(f"[{i}/{len(files)}] {f.name[:26]:<28} scanned={res.n_scanned:>7,} "
          f"triad={res.n_triad_complete:>3} written={res.n_written:>3}  {time.time()-t:.0f}s", flush=True)
print(f"TOTAL {time.time()-t0:.0f}s", flush=True)
with connect() as c:
    for r in c.execute("SELECT source_environment, COUNT(*) n FROM candidates GROUP BY 1 ORDER BY n DESC"):
        print(f"  {r[0]:<22} {r[1]}", flush=True)
