#!/usr/bin/env python3
"""Rebuild the reference structures for enzymes that just acquired a PDB deposit.

Run after `pipeline.structure.deposits.link()`. Those enzymes already have a structure on
disk, folded by ESMFold or fetched from AlphaFold; `reference.build(only=...)` reruns them
and `sources_for()` now finds a deposit at the head of the list, so the prediction is
replaced by the experimental coordinates and the row's source flips to 'pdb'.

Fetch-only for the ones that resolve, which is why this takes minutes rather than the
hours the bulk fold took: no ESMFold call is made for an enzyme whose deposit answers.
"""

from __future__ import annotations

import json
import sys

from pipeline.db import connect
from pipeline.structure import reference


def main() -> int:
    with connect() as c:
        rows = c.execute(
            "SELECT e.enzyme_id, r.source FROM characterised_enzymes e "
            "LEFT JOIN reference_structures r USING(enzyme_id) "
            "WHERE e.pdb_ids_json IS NOT NULL AND e.pdb_ids_json NOT IN ('', '[]') "
            "  AND (r.source IS NULL OR r.source != 'pdb') "
            "ORDER BY e.enzyme_id").fetchall()
    targets = [r[0] for r in rows]
    was = {r[0]: (r[1] or "none") for r in rows}
    if not targets:
        print("nothing to rebuild: every linked enzyme already uses its deposit")
        return 0

    print(f"rebuilding {len(targets)} enzymes now carrying a deposit", flush=True)
    report = reference.build(only=targets, label="deposit-relink")

    with connect() as c:
        after = {r[0]: r[1] for r in c.execute(
            f"SELECT enzyme_id, source FROM reference_structures "
            f"WHERE enzyme_id IN ({','.join('?' * len(targets))})", targets)}
    flipped = [e for e in targets if after.get(e) == "pdb" and was.get(e) != "pdb"]
    print(f"\n{len(flipped)} now use an experimental structure "
          f"(were: {json.dumps({s: sum(1 for e in flipped if was[e] == s) for s in set(was[e] for e in flipped)})})")
    for line in report["skipped"]:
        print(f"  SKIPPED  {line}")
    for line in report["failed"]:
        print(f"  FAILED   {line}")
    still = [e for e in targets if after.get(e) != "pdb"]
    if still:
        print(f"\n{len(still)} kept a prediction because no deposit returned coordinates:")
        for e in still[:20]:
            print(f"    {e}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
