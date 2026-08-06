#!/usr/bin/env python3
"""Lift the primary reference DOI out of free text into a column.

The PAZy import recorded each enzyme's citation as "Primary reference doi:10.xxxx/yyy"
inside `activity_substrate_notes`. Nothing was lost, but a DOI inside a sentence cannot be
grouped on, and grouping is the whole point: an ordinal ranking exists only where one paper
assayed several enzymes under one protocol, so "how many enzymes share a DOI" is the
question that decides whether that route is worth taking at all.
"""

from __future__ import annotations

import pathlib
import re
import sys
import urllib.parse
from collections import Counter

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline.db import connect

# DOIs end at whitespace or a sentence-final full stop; the trailing period in
# "doi:10.1021/acscatal.1c00126." is punctuation, not part of the identifier.
DOI_RE = re.compile(r"doi:\s*(10\.\d{4,9}/[^\s,;]+?)\.?(?:\s|$)", re.I)


def main() -> int:
    with connect() as c:
        rows = c.execute(
            "SELECT enzyme_id, activity_substrate_notes FROM characterised_enzymes "
            "WHERE activity_substrate_notes LIKE '%doi:%'").fetchall()
    updates, misses = [], []
    for enzyme_id, notes in rows:
        # One PAZy record stores the DOI percent-encoded ("10.1021%2Facs..."), which is a
        # URL fragment that escaped its URL. Unquoted first, or the slash never matches.
        text = urllib.parse.unquote(notes or "")
        m = DOI_RE.search(text)
        if m:
            updates.append((m.group(1).rstrip("."), enzyme_id))
        else:
            misses.append(enzyme_id)

    with connect() as c:
        c.executemany("UPDATE characterised_enzymes SET primary_doi=? WHERE enzyme_id=?",
                      updates)
        c.commit()
        counts = Counter(d for d, _ in updates)
        n_col = c.execute("SELECT COUNT(*) FROM characterised_enzymes "
                          "WHERE primary_doi IS NOT NULL").fetchone()[0]

    multi = {d: n for d, n in counts.items() if n > 1}
    in_multi = sum(multi.values())
    print(f"{len(rows)} rows carried a DOI in free text; {len(updates)} parsed, "
          f"{len(misses)} unparsed")
    print(f"{n_col} enzymes now carry primary_doi, across {len(counts)} distinct papers\n")
    print(f"papers covering MORE THAN ONE enzyme: {len(multi)}")
    print(f"enzymes inside such a paper:          {in_multi} "
          f"({100*in_multi/max(len(updates),1):.0f}% of those with a DOI)")
    print("\nlargest comparison sets, which are where a ranking can be read:")
    for doi, n in sorted(multi.items(), key=lambda kv: -kv[1])[:12]:
        print(f"   {n:3d} enzymes   {doi}")
    if misses:
        print(f"\nunparsed: {', '.join(misses[:6])}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
