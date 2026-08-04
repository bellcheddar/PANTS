"""UniProt REST client: the only route by which a real sequence enters PANTS.

Every sequence in `characterised_enzymes` is fetched from UniProt by accession or pulled
from a UniProt search. Nothing is ever typed in by hand: a mistyped residue in a seed
sequence would propagate silently into the HMM profiles, the embeddings and every score
downstream, and would be close to undetectable afterwards.

The REST API paginates through an RFC 5988 Link header rather than an offset parameter,
so `search()` follows `rel="next"` until exhausted.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from typing import Any, Dict, Iterator, List, Optional

from . import config, http

ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"

# Fields requested from search endpoints. Keep this tight: the default response carries
# the full cross-reference block, which is megabytes per entry at scale.
SEARCH_FIELDS = "accession,id,protein_name,organism_name,organism_id,length,sequence,lineage,reviewed"

_LINK_NEXT = re.compile(r'<([^>]+)>;\s*rel="next"')


@dataclass
class Entry:
    """One UniProt entry, flattened to what PANTS actually stores."""
    accession: str
    entry_name: Optional[str]
    protein_name: Optional[str]
    organism: Optional[str]
    taxid: Optional[int]
    sequence: str
    length: int
    lineage: Optional[str]
    reviewed: bool

    @property
    def is_plausible_protein(self) -> bool:
        """Reject anything with non-standard residues before it reaches an HMM or ESM-2."""
        return bool(self.sequence) and set(self.sequence) <= set("ACDEFGHIKLMNPQRSTVWY")


def _flatten(obj: Dict[str, Any]) -> Entry:
    seq = (obj.get("sequence") or {}).get("value", "")
    org = obj.get("organism") or {}
    protein = (((obj.get("proteinDescription") or {}).get("recommendedName") or {})
               .get("fullName") or {}).get("value")
    if not protein:
        subs = (obj.get("proteinDescription") or {}).get("submissionNames") or []
        if subs:
            protein = (subs[0].get("fullName") or {}).get("value")
    lineage = org.get("lineage")
    return Entry(
        accession=obj.get("primaryAccession", ""),
        entry_name=obj.get("uniProtkbId"),
        protein_name=protein,
        organism=org.get("scientificName"),
        taxid=org.get("taxonId"),
        sequence=seq,
        length=len(seq),
        lineage="; ".join(lineage) if isinstance(lineage, list) else lineage,
        reviewed=obj.get("entryType", "").startswith("UniProtKB reviewed"),
    )


def fetch(accession: str) -> Optional[Entry]:
    """One entry by accession. Returns None if UniProt does not know it (404)."""
    obj = http.get_json(ENTRY_URL.format(accession=accession))
    return _flatten(obj) if obj else None


def fetch_many(accessions: List[str]) -> Dict[str, Optional[Entry]]:
    """Sequential fetch. Deliberately not parallel: this runs over ~10^2 seed accessions
    once, and hammering UniProt to save twenty seconds is not a good trade."""
    return {acc: fetch(acc) for acc in accessions}


def search(query: str, size: int = 500, max_results: Optional[int] = None,
           fields: str = SEARCH_FIELDS) -> Iterator[Entry]:
    """Yield entries for a UniProt query string, following pagination to exhaustion.

    `max_results` caps the walk: the ESTHER negative families run to tens of thousands of
    entries and the harvest samples rather than taking everything.
    """
    url = config.UNIPROT_REST_URL
    params: Optional[Dict[str, Any]] = {
        "query": query, "format": "json", "size": min(size, 500), "fields": fields,
    }
    seen = 0
    while url:
        resp = http.get(url, params=params)
        resp.raise_for_status()
        for obj in resp.json().get("results", []):
            yield _flatten(obj)
            seen += 1
            if max_results is not None and seen >= max_results:
                return
        # Subsequent pages come fully-formed from the Link header: passing params again
        # would clobber the cursor and restart the walk from page one, forever.
        m = _LINK_NEXT.search(resp.headers.get("Link", ""))
        url, params = (m.group(1) if m else None), None


def count(query: str) -> Optional[int]:
    """Total hits for a query, from the x-total-results header, without downloading them."""
    resp = http.get(config.UNIPROT_REST_URL,
                    params={"query": query, "format": "json", "size": 0})
    resp.raise_for_status()
    total = resp.headers.get("x-total-results")
    return int(total) if total is not None else None
