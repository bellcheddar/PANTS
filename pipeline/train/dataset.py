"""Assemble the training set, and mark what is excluded from it and why.

Exclusions are recorded on the row rather than applied by deleting it. The catalogue is a
deliverable in its own right (spec section 7), so a fragment still belongs in it; what it
must not do is quietly become a training example. Marking keeps the decision auditable:
`SELECT exclusion_reason, COUNT(*) ... GROUP BY 1` answers "what did you throw away?"
without re-running anything.

Two filters, both anchored to evidence rather than intuition:

  Fragments come from UniProt's own `fragment:true` flag, not a length threshold. A short
  sequence is not necessarily a fragment and a fragment is not necessarily short.

  That flag is now asked for EVERY accession in the catalogue. It used to be asked only of
  `ec:3.1.1.101`, which is the PET hydrolase EC number, so a fragment arriving through PAZy
  or ESTHER under any other EC was invisible unless it also happened to fall outside the
  length window -- caught by the right rule for the wrong reason, and not caught at all if
  it was fragmentary but of ordinary length.

  The length window is derived from the EXPERIMENTALLY EVIDENCED positives (262 to 319 aa
  for the PET hydrolases, excluding MHETase which is a different fold entirely), padded
  generously. The EC-auto-annotated set stretches 63 to 835 aa, which drags in fragments
  and probable multi-domain proteins that would otherwise teach the model that "polyesterase"
  spans an order of magnitude in size.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Set, Tuple

from .. import http, uniprot
from ..db import connect, now

# Padded around the evidenced range (262 to 319 aa). Wide enough to keep genuine
# variation, narrow enough to exclude fragments and multi-domain proteins.
LENGTH_MIN = 200
LENGTH_MAX = 450

FRAGMENT_QUERY = "ec:3.1.1.101 AND fragment:true"

# UniProt rejects a query string past a few kB, so accessions go up in batches.
_ACCESSION_BATCH = 80


def catalogue_accessions() -> List[str]:
    with connect() as conn:
        return [r[0].split("-")[0] for r in conn.execute(
            "SELECT DISTINCT uniprot FROM characterised_enzymes "
            "WHERE uniprot IS NOT NULL AND uniprot != ''")]


def fetch_fragment_accessions(accessions: Optional[List[str]] = None) -> Set[str]:
    """Accessions UniProt itself flags as fragments.

    Asks about the accessions actually in the catalogue rather than about an EC number.
    `fragment:true` is UniProt's own filter, so this stays their definition and not a
    length heuristic wearing their name; combining it with the accession list just narrows
    who is asked about.
    """
    accs = accessions if accessions is not None else catalogue_accessions()
    if not accs:
        return set()
    found: Set[str] = set()
    for i in range(0, len(accs), _ACCESSION_BATCH):
        chunk = accs[i:i + _ACCESSION_BATCH]
        q = "(" + " OR ".join(f"accession:{a}" for a in chunk) + ") AND fragment:true"
        found.update(e.accession for e in uniprot.search(q, max_results=len(chunk)))
    return found


def apply_filters(length_min: int = LENGTH_MIN, length_max: int = LENGTH_MAX
                  ) -> Dict[str, int]:
    """Mark rows excluded from training. Idempotent: recomputes every flag from scratch."""
    fragments = fetch_fragment_accessions()
    counts: Dict[str, int] = {}

    with connect() as conn:
        conn.execute(
            "UPDATE characterised_enzymes SET excluded_from_training=0, exclusion_reason=NULL")

        # UniProt's own fragment flag.
        if fragments:
            marks = ",".join("?" * len(fragments))
            conn.execute(
                f"UPDATE characterised_enzymes SET is_fragment=1, excluded_from_training=1, "
                f"exclusion_reason='UniProt Fragment flag' WHERE uniprot IN ({marks})",
                tuple(fragments))

        # Length window, applied only to sequences that have one.
        conn.execute(
            "UPDATE characterised_enzymes SET excluded_from_training=1, "
            "exclusion_reason='length outside " f"{length_min}-{length_max}" " aa' "
            "WHERE excluded_from_training=0 AND seq_length IS NOT NULL "
            "AND (seq_length < ? OR seq_length > ?)", (length_min, length_max))

        # No sequence at all: the variants whose mutation sets were never confirmed.
        conn.execute(
            "UPDATE characterised_enzymes SET excluded_from_training=1, "
            "exclusion_reason='no sequence (mutation set unconfirmed)' "
            "WHERE excluded_from_training=0 AND sequence IS NULL")

        for r in conn.execute(
                "SELECT COALESCE(exclusion_reason,'included') reason, COUNT(*) n "
                "FROM characterised_enzymes GROUP BY 1 ORDER BY n DESC"):
            counts[r["reason"]] = r["n"]
    return counts


def training_set(include_near_misses: bool = True) -> Tuple[List[Tuple[str, str]], List[int], List[str]]:
    """(records, labels, evidence_tiers) for the head.

    Near misses are labelled NEGATIVE: they are triad-bearing esterases active on soluble
    esters but not meaningful degraders of crystalline PET, which is precisely the boundary
    the head has to learn (spec section 5.2). Including them as negatives is the whole
    point of having curated them.
    """
    where_nm = "" if include_near_misses else " AND is_near_miss=0"
    with connect() as conn:
        rows = list(conn.execute(
            "SELECT enzyme_id, sequence, is_positive, is_negative, is_near_miss, source_ref "
            "FROM characterised_enzymes "
            "WHERE excluded_from_training=0 AND sequence IS NOT NULL "
            "AND family != 'mhetase_like'" + where_nm))
    records = [(r["enzyme_id"], r["sequence"]) for r in rows]
    labels = [1 if r["is_positive"] else 0 for r in rows]
    tiers = [r["source_ref"] or "unknown" for r in rows]
    return records, labels, tiers
