#!/usr/bin/env python3
"""Populate the date UniProt first made each catalogue entry public.

The evaluation protocol calls for a prospective holdout: train on what was known before a
cutoff, test on what appeared after it, which is the only split that answers "would this
have found the enzymes we now know about". `pdb_release_date` cannot serve, because it is
empty for the entire catalogue -- most of these enzymes have never been crystallised.

UniProt's first-public date exists for anything with an accession and is the earliest
defensible proxy for when the sequence entered the public record. It is a proxy and not a
discovery date: an enzyme can be characterised years before its sequence is deposited, and
a 2016 accession for IsPETase does not mean nobody had it in 2015. What it does give is a
monotone ordering that was not derived from the labels, which is what the split needs.
"""

from __future__ import annotations

import json
import pathlib
import sys
import urllib.parse
from typing import Dict, List

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import http
from pipeline.db import connect

SEARCH = "https://rest.uniprot.org/uniprotkb/search"
BATCH = 80


def fetch(accessions: List[str]) -> Dict[str, str]:
    out: Dict[str, str] = {}
    for i in range(0, len(accessions), BATCH):
        chunk = accessions[i:i + BATCH]
        q = " OR ".join(f"accession:{a}" for a in chunk)
        data = http.get_json(f"{SEARCH}?format=json&size=500&fields=accession,date_created"
                             f"&query={urllib.parse.quote(q)}")
        for r in (data or {}).get("results", []):
            d = (r.get("entryAudit") or {}).get("firstPublicDate")
            if d:
                out[r["primaryAccession"]] = d
        print(f"  {min(i + BATCH, len(accessions))}/{len(accessions)} scanned, "
              f"{len(out)} dated", flush=True)
    return out


def main() -> int:
    with connect() as c:
        accs = sorted({r[0].split("-")[0] for r in c.execute(
            "SELECT DISTINCT uniprot FROM characterised_enzymes "
            "WHERE uniprot IS NOT NULL AND uniprot != ''")})
    dates = fetch(accs)
    with connect() as c:
        c.executemany(
            "UPDATE characterised_enzymes SET uniprot_first_public=? "
            "WHERE uniprot LIKE ? || '%'",
            [(d, a) for a, d in dates.items()])
        c.commit()
        n = c.execute("SELECT COUNT(*) FROM characterised_enzymes "
                      "WHERE uniprot_first_public IS NOT NULL").fetchone()[0]
        span = c.execute("SELECT MIN(uniprot_first_public), MAX(uniprot_first_public) "
                         "FROM characterised_enzymes "
                         "WHERE uniprot_first_public IS NOT NULL").fetchone()
    print(f"\n{len(dates)} of {len(accs)} accessions dated; {n} catalogue rows carry a date")
    print(f"spanning {span[0]} to {span[1]}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
