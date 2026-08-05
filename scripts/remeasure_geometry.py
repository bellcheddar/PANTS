"""Re-measure active-site geometry for every structure already on disk.

Needed whenever a measurement changes. Folding is the expensive step and the coordinates
do not change, so re-measuring is cheap and must never require re-folding.
"""

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

import json

from pipeline import config
from pipeline.db import connect, now, retry_write
from pipeline.structure import geometry

STRUCT_DIR = config.APP_DIR / "static" / "structures"

if __name__ == "__main__":
    with connect() as c:
        rows = list(c.execute(
            "SELECT s.candidate_id, s.mmcif_path FROM structures s ORDER BY s.candidate_id"))
    print(f"{len(rows)} structures to re-measure", flush=True)

    changed = missing = failed = 0
    for i, (cid, name) in enumerate(rows, 1):
        path = STRUCT_DIR / name
        if not path.exists():
            for alt in (STRUCT_DIR / f"{cid}.pdb", STRUCT_DIR / f"{cid}.cif"):
                if alt.exists():
                    path = alt
                    break
            else:
                missing += 1
                continue
        try:
            site = geometry.measure(path)
        except Exception as exc:
            failed += 1
            print(f"  {cid}: {type(exc).__name__} {exc}", flush=True)
            continue

        def _do(cid=cid, site=site):
            with connect() as c:
                c.execute(
                    "UPDATE geometry SET oxyanion_n1_dist_A=?, oxyanion_n2_dist_A=?, "
                    " oxyanion_n1_resnum=?, oxyanion_n2_resnum=?, oxyanion_n2_angle_deg=? "
                    "WHERE candidate_id=?",
                    (site.oxyanion_n1_A, site.oxyanion_n2_A, site.oxyanion_n1_resnum,
                     site.oxyanion_n2_resnum, site.oxyanion_n2_angle_deg, cid))
        retry_write(_do)
        changed += 1
        if i % 50 == 0:
            print(f"  {i}/{len(rows)}", flush=True)

    print(f"\nre-measured {changed}, missing file {missing}, failed {failed}")
    with connect() as c:
        n = c.execute("SELECT COUNT(*) FROM geometry WHERE oxyanion_n2_resnum IS NOT NULL").fetchone()[0]
        tot = c.execute("SELECT COUNT(*) FROM geometry").fetchone()[0]
        print(f"rows with an identified second donor: {n}/{tot}")
