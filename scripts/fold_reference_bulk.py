"""Build structures for the bulk characterised set, in resumable batches.

The named lineages have structures; the ~340 PAZy-derived enzymes the head is actually
trained on do not, so they cannot appear in any structural view. This closes that.

Resumable by construction: anything already in reference_structures is skipped, so an
interrupted run continues where it stopped rather than starting over. That matters because
the work is hours long and the machine it runs on is also somebody's laptop.

AlphaFold is tried before ESMFold wherever a UniProt accession exists, which is most of
them -- a download is seconds against roughly ninety for a fold, and a model built by
AlphaFold from an MSA is generally the better structure anyway.
"""

from __future__ import annotations

import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline.db import connect
from pipeline.structure import fold, reference

BATCH = 25


def pending() -> list:
    with connect() as c:
        return [r[0] for r in c.execute(
            "SELECT enzyme_id FROM characterised_enzymes ce "
            "WHERE ce.enzyme_id LIKE 'PAZy:%' AND ce.sequence IS NOT NULL "
            "  AND ce.seq_length <= ? "
            "  AND NOT EXISTS (SELECT 1 FROM reference_structures rs "
            "                  WHERE rs.enzyme_id = ce.enzyme_id) "
            "ORDER BY ce.enzyme_id", (fold.MAX_FOLD_LENGTH,))]


if __name__ == "__main__":
    todo = pending()
    print(f"{len(todo)} structures to build", flush=True)
    t0 = time.time()
    built = failed = 0
    while todo:
        chunk, todo = todo[:BATCH], todo[BATCH:]
        rep = reference.build(only=chunk, label="pazy-bulk")
        built += len(rep["built"])
        failed += len(rep["failed"])
        for f in rep["failed"]:
            print(f"  FAILED {f[:110]}", flush=True)
        for s in rep["skipped"]:
            if "KNOCKOUT" in s:
                print(f"  {s[:130]}", flush=True)
        el = time.time() - t0
        rate = built / el if el else 0
        print(f"  {built} built, {failed} failed, {len(todo)} left  "
              f"({el/60:.0f} min elapsed"
              + (f", ~{len(todo)/rate/60:.0f} min remaining)" if rate else ")"), flush=True)
    print(f"\nDONE: {built} built, {failed} failed in {(time.time()-t0)/60:.0f} min")
