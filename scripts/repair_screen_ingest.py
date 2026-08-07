#!/usr/bin/env python3
"""Repair two things the screen ingests did to rows that already existed.

BOTH were caused by the same `ON CONFLICT DO UPDATE` clause treating an existing row as
something to overwrite rather than something to add to.

1. PROVENANCE DESTROYED. `activity_substrate_notes=excluded.activity_substrate_notes`
   replaced the note wholesale for the 189 enzymes the screens matched by sequence, so an
   entry that recorded "PAZy 409: measured activity on PET, primary reference doi:..." now
   records only the screen result. The evidence was not contradicted, it was deleted.
   Restored from the published v0.2.0 Zenodo deposit -- which is what a versioned deposit
   is for -- with the screen note APPENDED rather than substituted.

2. CONTRADICTORY LABELS. Nine enzymes ended up flagged `is_positive=1` AND
   `within_family_basis='measured-inactive'`, and five of them were being drawn into both
   the positive and the negative side of the same cross-validation. Two distinct cases, and
   they do not get the same treatment:

   a. Three were positive only by EC number, assigned automatically by similarity
      (ECO:0000256), and a screen then measured 0 +/- 0 uM of product. A measurement beats
      an annotation; that is the premise the whole project rests on. They become negatives.

   b. Six are PAZy entries recording MEASURED PET activity -- and four of those cite the
      very paper whose data classed them inactive here. That is not two experiments
      disagreeing. It is my threshold being stricter than the curator's reading of one
      experiment: 0.090% against a 0.1% floor, or 288 +/- 311 uM where the standard
      deviation exceeds the value. Manufacturing a negative out of a threshold disagreement
      and then testing a model against it would be inventing evidence. Their
      measured-inactive basis is withdrawn, they are marked CONTESTED and excluded from
      training, and the borderline number is kept on the row so the call can be revisited.
"""

from __future__ import annotations

import csv
import io
import pathlib
import sys
import urllib.request
from typing import Dict

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline.db import connect

SNAPSHOT = ("https://zenodo.org/api/records/21827812/files/"
            "characterised_enzymes.csv/content")
SCREEN_MARKERS = ("ACS Catalysis 2025 high-throughput screen",
                  "Science 2025 PET depolymerase landscape")


def fetch_snapshot() -> Dict[str, dict]:
    with urllib.request.urlopen(SNAPSHOT, timeout=120) as fh:
        text = fh.read().decode("utf-8", "replace")
    return {r["enzyme_id"]: r for r in csv.DictReader(io.StringIO(text))}


def restore_notes(snapshot: Dict[str, dict]) -> int:
    with connect() as c:
        rows = c.execute(
            "SELECT enzyme_id, activity_substrate_notes, source_ref "
            "FROM characterised_enzymes "
            "WHERE source_ref NOT IN ('ACS-screen-measured','Science-landscape-measured')"
        ).fetchall()
        updates = []
        for eid, notes, _ in rows:
            if not notes or not any(m in notes for m in SCREEN_MARKERS):
                continue
            original = (snapshot.get(eid) or {}).get("activity_substrate_notes", "").strip()
            if not original or original in notes:
                continue
            updates.append((f"{original} SUBSEQUENTLY SCREENED: {notes}", eid))
        c.executemany("UPDATE characterised_enzymes SET activity_substrate_notes=? "
                      "WHERE enzyme_id=?", updates)
        c.commit()
    return len(updates)


def resolve_contradictions() -> dict:
    report = {"annotation_overruled": [], "contested": []}
    with connect() as c:
        rows = c.execute(
            "SELECT enzyme_id, source_ref, activity_substrate_notes "
            "FROM characterised_enzymes "
            "WHERE is_positive=1 AND within_family_basis='measured-inactive'").fetchall()
        for eid, src, notes in rows:
            annotation_only = src == "EC-auto-annotated"
            if annotation_only:
                # A measurement beats a similarity annotation.
                c.execute(
                    "UPDATE characterised_enzymes SET is_positive=0, is_near_miss=1, "
                    " activity_substrate_notes = activity_substrate_notes || "
                    " ' RESOLUTION: the EC 3.1.1.101 label was assigned automatically by "
                    "similarity and a screen then measured no product release, so the "
                    "measurement governs and this is recorded as a within-family negative.' "
                    "WHERE enzyme_id=?", (eid,))
                report["annotation_overruled"].append(eid)
            else:
                # PAZy recorded a published, MEASURED PET activity. Withdraw the negative
                # rather than assert it over the source.
                c.execute(
                    "UPDATE characterised_enzymes SET within_family_basis=NULL, "
                    " is_near_miss=0, excluded_from_training=1, "
                    " exclusion_reason='contested: published PET activity vs a screen "
                    "result below this threshold', "
                    " activity_substrate_notes = activity_substrate_notes || "
                    " ' RESOLUTION: PAZy records measured, published PET activity for this "
                    "enzyme while a screen result falls below the activity threshold used "
                    "here -- in several cases the SAME paper. That is a threshold "
                    "disagreement rather than a contradicting experiment, so no negative is "
                    "claimed and the enzyme is held out of training as contested.' "
                    "WHERE enzyme_id=?", (eid,))
                report["contested"].append(eid)
        c.commit()
    return report


def main() -> int:
    print("fetching the published v0.2.0 snapshot ...", flush=True)
    snapshot = fetch_snapshot()
    print(f"  {len(snapshot)} rows\n")

    n = restore_notes(snapshot)
    print(f"provenance restored on {n} rows (original note kept, screen note appended)\n")

    r = resolve_contradictions()
    print(f"annotation overruled by measurement : {len(r['annotation_overruled'])}")
    for e in r["annotation_overruled"]:
        print(f"    {e}")
    print(f"contested, no negative claimed      : {len(r['contested'])}")
    for e in r["contested"]:
        print(f"    {e}")

    with connect() as c:
        left = c.execute("SELECT COUNT(*) FROM characterised_enzymes "
                         "WHERE is_positive=1 AND within_family_basis='measured-inactive'"
                         ).fetchone()[0]
        neg = c.execute("SELECT COUNT(*) FROM characterised_enzymes "
                        "WHERE within_family_basis='measured-inactive'").fetchone()[0]
    print(f"\nrows still on both sides: {left}   (must be 0)")
    print(f"measured-inactive negatives now: {neg}")
    return 0 if left == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
