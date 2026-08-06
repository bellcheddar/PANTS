#!/usr/bin/env python3
"""Give a bulk-imported enzyme its published name, everywhere at once.

Micpa-PETase and Kutbu-PETase arrived through the PAZy import as `PAZy:270` and `PAZy:276`
with PAZy's own abbreviations, `Mipa-P` and `Kubu-P`. They are not bulk entries in any
meaningful sense -- both are named in the Science landscape paper, both have crystal
structures (8YTU at 1.34 A and 8YTW at 2.65 A), and both have a measured triad. What kept
them out of the reference tables, the lineage pages and the enzyme deep-dives was the
identifier they happened to be imported under: every one of those views filters on
`enzyme_id NOT LIKE 'PAZy:%'`.

So this renames the row rather than duplicating it. A second row carrying the same
sequence under a nicer name would double-count the enzyme in every total on the site, and
the totals are the thing this project asks to be trusted on.

The rename has to cascade by hand. `reference_structures` and `activity_measurements`
reference `characterised_enzymes(enzyme_id)`, and `reference_geometry` references
`reference_structures`, all with `ON UPDATE NO ACTION` -- so SQLite will not carry the
change down, and with `foreign_keys=ON` it refuses the parent update outright while
children point at the old value. Children are therefore updated first, inside one
transaction, with the constraint deferred for its duration.

The coordinate file is renamed alongside, because `coord_path` is what the viewer fetches
and a stale path is a 404 in a 3D panel rather than an error anyone would see.
"""

from __future__ import annotations

import argparse
import pathlib
import sys
from typing import List, Tuple

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config
from pipeline.db import connect

REF_DIR = config.STATIC_DIR / "reference_structures"

# (old id, new id, common name, one-line headline for the summary table)
PROMOTIONS: List[Tuple[str, str, str, str]] = [
    ("PAZy:270", "Micpa-PETase", "Micpa-PETase",
     "Actinomycete PETase from Micromonospora pattaloongensis, solved at 1.34 A"),
    ("PAZy:276", "Kutbu-PETase", "Kutbu-PETase",
     "Actinomycete PETase from Kutzneria buriramensis, solved at 2.65 A"),
]


def _safe(enzyme_id: str) -> str:
    """The filename rule reference.build uses, so the renamed file is the one it expects."""
    return enzyme_id.replace("/", "_").replace("*", "s")


def promote(old: str, new: str, common: str, headline: str, dry_run: bool = False) -> bool:
    with connect() as c:
        row = c.execute("SELECT enzyme_id FROM characterised_enzymes WHERE enzyme_id=?",
                        (old,)).fetchone()
        if not row:
            existing = c.execute("SELECT enzyme_id FROM characterised_enzymes WHERE enzyme_id=?",
                                 (new,)).fetchone()
            print(f"  {old} -> {new}: {'already promoted' if existing else 'NOT FOUND'}")
            return bool(existing)
        if c.execute("SELECT 1 FROM characterised_enzymes WHERE enzyme_id=?", (new,)).fetchone():
            print(f"  {old} -> {new}: target id already taken, refusing")
            return False

        coord = c.execute("SELECT coord_path FROM reference_structures WHERE enzyme_id=?",
                          (old,)).fetchone()
        old_file = REF_DIR / coord[0] if coord and coord[0] else None
        new_name = f"{_safe(new)}.pdb"

        if dry_run:
            print(f"  {old} -> {new}  (would rename {coord[0] if coord else 'no structure'} "
                  f"to {new_name})")
            return True

        c.execute("PRAGMA defer_foreign_keys=ON")
        # Children first: with foreign_keys ON, updating the parent while children still
        # point at the old id is the case SQLite rejects.
        c.execute("UPDATE reference_geometry   SET enzyme_id=? WHERE enzyme_id=?", (new, old))
        c.execute("UPDATE reference_structures SET enzyme_id=?, coord_path=? WHERE enzyme_id=?",
                  (new, new_name, old))
        c.execute("UPDATE activity_measurements SET enzyme_id=? WHERE enzyme_id=?", (new, old))
        c.execute("UPDATE characterised_enzymes SET enzyme_id=?, common_name=?, headline=? "
                  "WHERE enzyme_id=?", (new, common, headline, old))
        # Anything naming it in prose: the notes still read "PAZy Mipa-P: ...".
        c.execute("UPDATE characterised_enzymes "
                  "SET activity_substrate_notes = REPLACE(activity_substrate_notes, ?, ?) "
                  "WHERE enzyme_id=?", (f"PAZy {common.split('-')[0][:4]}", f"PAZy {common}", new))
        c.commit()

    if old_file and old_file.exists():
        old_file.rename(REF_DIR / new_name)
    print(f"  {old} -> {new}   coordinates now {new_name}")
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    ok = 0
    for old, new, common, headline in PROMOTIONS:
        ok += bool(promote(old, new, common, headline, args.dry_run))

    if not args.dry_run:
        with connect() as c:
            for _, new, _, _ in PROMOTIONS:
                r = c.execute(
                    "SELECT ce.enzyme_id, ce.common_name, ce.uniprot, ce.organism, "
                    "       rs.source_id, rs.resolution_A, rg.triad_ser_resnum "
                    "FROM characterised_enzymes ce "
                    "LEFT JOIN reference_structures rs ON rs.enzyme_id=ce.enzyme_id "
                    "LEFT JOIN reference_geometry rg ON rg.enzyme_id=ce.enzyme_id "
                    "WHERE ce.enzyme_id=?", (new,)).fetchone()
                if r:
                    print(f"    {r[0]}: {r[3]}, PDB {r[4]} at {r[5]} A, catalytic Ser{r[6]}")
    print(f"\n{ok} of {len(PROMOTIONS)} promoted")
    return 0


if __name__ == "__main__":
    sys.exit(main())
