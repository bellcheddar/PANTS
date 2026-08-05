"""Keep folding new candidates as recall produces them.

The first fold run snapshotted its candidate list at start, so everything recall found
afterwards was never queued. This drains continuously instead: one model load, then a
loop that folds whatever has no structure yet.

Waits for any existing fold process to exit first. Two ESMFold processes would hold two
copies of an 8.4 GB model (~11 GB resident each) and compete for the same cores.
"""
import pathlib
import sys, time, subprocess
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from pathlib import Path
from pipeline import config
from pipeline.db import connect
from pipeline.db.manifest import stage_manifest
from pipeline.structure import fold

def pending():
    with connect() as c:
        return [(r[0], r[1]) for r in c.execute(
            "SELECT c.candidate_id, c.sequence FROM candidates c "
            "LEFT JOIN structures s ON s.candidate_id = c.candidate_id "
            "WHERE s.candidate_id IS NULL AND COALESCE(c.structure_deferred,0)=0 "
            "AND c.seq_length <= ? ORDER BY c.recall_bitscore DESC",
            (fold.MAX_FOLD_LENGTH,))]

def other_fold_running():
    """True only if a PYTHON process is running run_fold.py.

    `pgrep -f run_fold.py` is not that test: it matches ANY command line containing the
    string, including the monitor shell that watches the job. That is exactly what
    happened here, and the drainer waited on an already-finished job while 78 candidates
    sat queued. Match on the executable being python, not on the string appearing.
    """
    out = subprocess.run(["ps", "-eo", "comm,args"], capture_output=True, text=True)
    for line in out.stdout.splitlines():
        parts = line.split(None, 1)
        if len(parts) == 2 and "python" in parts[0].lower() and "run_fold.py" in parts[1]:
            return True
    return False

while other_fold_running():
    print("  waiting for the existing fold run to finish...", flush=True)
    time.sleep(120)

print("loading ESMFold once for the whole drain", flush=True)
t0 = time.time()
tok, model = fold.load_model()
print(f"  model loaded in {time.time()-t0:.0f}s", flush=True)

GUT_DONE = config.INTERIM_DIR / "gut_done.txt"
idle_rounds = 0
total = 0
while True:
    todo = pending()
    if todo:
        idle_rounds = 0
        print(f"[drain] {len(todo)} pending", flush=True)
        for cid, seq in todo:
            t = time.time()
            try:
                pdb, plddt = fold.fold_one(tok, model, seq)
            except Exception as exc:
                print(f"  {cid} FAILED: {type(exc).__name__}: {exc}", flush=True)
                continue
            with fold._TMALIGN_LOCK:
                cif, rmsd, frac = fold.superpose_onto_reference(
                    pdb, config.INTERIM_DIR / "pdb" / "6EQE.cif")
            out = config.STRUCTURE_DIR
            path = out / f"{cid}.cif"
            if cif:
                path.write_text(cif)
                fold.write_viewer_pdb(cif, out / f"{cid}.pdb")
            else:
                path = out / f"{cid}.pdb"; path.write_text(pdb)
            site = fold.geometry.measure(path)
            fold._persist(cid, path, plddt, frac, rmsd, site)
            total += 1
            print(f"  {cid} {len(seq)}aa {time.time()-t:.0f}s pLDDT {plddt:.1f} "
                  f"triad {'yes' if site.triad_is_connected else 'NO'} "
                  f"cleft {site.cleft_width_A}", flush=True)
    else:
        idle_rounds += 1
        scan_done = GUT_DONE.exists() and len(GUT_DONE.read_text().split()) >= 50
        # Only stop once the scan that feeds this has finished AND nothing new has
        # appeared for a while: exiting on an empty queue alone would quit during a
        # gap between recall runs.
        if scan_done and idle_rounds >= 3:
            print(f"DRAIN COMPLETE: {total} structures added, recall finished, queue empty", flush=True)
            break
        print(f"  queue empty (scan_done={scan_done}), waiting...", flush=True)
        time.sleep(300)
