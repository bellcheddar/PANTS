#!/usr/bin/env python3
"""Which of the multi-enzyme papers can actually be read, and which need to be supplied.

An ordinal ranking exists only where one protocol assayed several enzymes, so the papers
worth reading are exactly those covering more than one enzyme in the catalogue: 31 of them,
holding 285 enzymes between them. Before extracting anything, this establishes what is
reachable, because a plan that assumes access to a paywalled Science supplement is a plan
that stops halfway through with no warning.

Three sources, cheapest first and none of them scraping a publisher:

  Europe PMC   full text as XML for anything in the open-access subset, which is the only
               route that gives machine-readable tables rather than a PDF to squint at.
  Unpaywall    the canonical answer to "is there a legal free copy, and where", including
               author manuscripts in repositories the publisher does not host.
  Crossref     always answers, and gives the title, journal and year needed to hand a
               human a list of exactly what to fetch.

Reports per paper: how many enzymes it covers, whether the full text is machine-readable,
and if not, where a copy is. Nothing is downloaded here beyond metadata; this is the survey
that decides where the effort goes.
"""

from __future__ import annotations

import json
import pathlib
import sys
import time
import urllib.parse
from typing import Any, Dict, List, Optional

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config, http
from pipeline.db import connect

EPMC_SEARCH = "https://www.ebi.ac.uk/europepmc/webservices/rest/search"
UNPAYWALL = "https://api.unpaywall.org/v2/{doi}?email={email}"
CROSSREF = "https://api.crossref.org/works/{doi}"
EMAIL = "marc@marcdeller.com"

OUT = config.INTERIM_DIR / "paper_access.json"


def multi_enzyme_papers() -> List[Dict[str, Any]]:
    """Papers covering more than one catalogued enzyme, largest panel first."""
    with connect() as c:
        rows = c.execute(
            "SELECT primary_doi, COUNT(*) n, "
            "       SUM(CASE WHEN is_positive=1 THEN 1 ELSE 0 END) n_pos, "
            "       SUM(CASE WHEN within_family_basis IS NOT NULL THEN 1 ELSE 0 END) n_wfn "
            "FROM characterised_enzymes WHERE primary_doi IS NOT NULL "
            "GROUP BY 1 HAVING n > 1 ORDER BY n DESC").fetchall()
    return [{"doi": r[0], "n_enzymes": r[1], "n_positive": r[2],
             "n_within_family_negative": r[3]} for r in rows]


def europepmc(doi: str) -> Dict[str, Any]:
    """Identifiers and, crucially, whether the machine-readable full text exists."""
    q = urllib.parse.quote(f'DOI:"{doi}"')
    data = http.get_json(f"{EPMC_SEARCH}?query={q}&format=json&resultType=core&pageSize=1")
    hits = ((data or {}).get("resultList") or {}).get("result") or []
    if not hits:
        return {}
    r = hits[0]
    return {"pmid": r.get("pmid"), "pmcid": r.get("pmcid"),
            "is_open_access": r.get("isOpenAccess") == "Y",
            # hasTextMinedTerms is not the same thing: the XML is what carries the tables.
            "has_full_text_xml": r.get("hasTextMinedTerms") == "Y" or bool(r.get("pmcid")),
            "has_supplementary": r.get("hasSuppl") == "Y",
            "title": r.get("title"), "journal": (r.get("journalInfo") or {}).get("journal", {}).get("title"),
            "year": r.get("pubYear")}


def unpaywall(doi: str) -> Dict[str, Any]:
    data = http.get_json(UNPAYWALL.format(doi=urllib.parse.quote(doi), email=EMAIL))
    if not data:
        return {}
    loc = data.get("best_oa_location") or {}
    return {"is_oa": data.get("is_oa"), "oa_status": data.get("oa_status"),
            "oa_url": loc.get("url_for_pdf") or loc.get("url"),
            "host_type": loc.get("host_type")}


def crossref(doi: str) -> Dict[str, Any]:
    data = http.get_json(CROSSREF.format(doi=urllib.parse.quote(doi)))
    m = (data or {}).get("message") or {}
    return {"title": (m.get("title") or [None])[0],
            "journal": (m.get("container-title") or [None])[0],
            "year": (m.get("issued", {}).get("date-parts") or [[None]])[0][0]}


def main() -> int:
    papers = multi_enzyme_papers()
    total = sum(p["n_enzymes"] for p in papers)
    print(f"{len(papers)} papers cover more than one catalogued enzyme, "
          f"{total} enzymes between them\n")

    for p in papers:
        doi = p["doi"]
        p.update(crossref(doi))
        e = europepmc(doi)
        p.update({k: v for k, v in e.items() if v is not None and k not in ("title",)})
        if not e.get("pmcid"):
            p.update(unpaywall(doi))
        # Machine-readable means the tables can be parsed; a PDF behind a login cannot.
        p["reachable"] = bool(e.get("pmcid")) or bool(p.get("oa_url"))
        p["machine_readable"] = bool(e.get("pmcid"))
        time.sleep(0.2)          # three public APIs, none of which owes us a firehose

    reach = [p for p in papers if p["reachable"]]
    xml = [p for p in papers if p["machine_readable"]]
    gap = [p for p in papers if not p["reachable"]]
    n = lambda rows: sum(r["n_enzymes"] for r in rows)

    print(f"{'enz':>4}  {'src':<14} {'year':<5} paper")
    for p in papers:
        src = ("EuropePMC XML" if p["machine_readable"]
               else (p.get("oa_status") or "OA") + " pdf" if p["reachable"] else "NOT REACHABLE")
        print(f"{p['n_enzymes']:>4}  {src:<14} {str(p.get('year') or '?'):<5} "
              f"{(p.get('title') or p['doi'])[:78]}")

    print(f"\nmachine-readable full text : {len(xml):>2} papers, {n(xml):>3} enzymes")
    print(f"reachable but as PDF only  : {len(reach)-len(xml):>2} papers, {n(reach)-n(xml):>3} enzymes")
    print(f"NOT reachable openly       : {len(gap):>2} papers, {n(gap):>3} enzymes")
    if gap:
        print("\n  to be supplied -- these are the ones worth your institutional access, "
              "largest panel first:")
        for p in gap:
            print(f"    {p['n_enzymes']:>3} enzymes  {p['doi']}")
            print(f"                 {(p.get('title') or '')[:90]}")
            print(f"                 {p.get('journal') or ''} {p.get('year') or ''}")

    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(json.dumps(papers, indent=2))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
