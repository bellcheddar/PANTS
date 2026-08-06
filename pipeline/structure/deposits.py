"""Find the experimental deposit for an enzyme that only had an accession.

The reference builder has always preferred a crystal structure over a prediction --
`sources_for()` puts every entry of `pdb_ids_json` ahead of AlphaFold. The preference
simply never fired for most of the reference set, because the PAZy import records a
UniProt accession and no PDB identifiers at all, so `pdb_ids_json` was empty for all 312
of its enzymes and every one of them fell through to a model. LCC, Cut190, TfCut2,
Est119, Thc_Cut1 and two dozen others were being shown an AlphaFold prediction while
their published crystal structures sat one cross-reference away.

This module closes that gap: it reads UniProt's PDB cross-references for an accession,
then ranks the deposits so the one written into `pdb_ids_json` is the one that actually
represents the enzyme.

Ranking is the whole problem. An accession maps to everything ever deposited for that
protein -- IsPETase alone has 55 entries, most of them engineered variants, several of
them deliberately inactivated. Taking the first, or the best resolution, gets you a
mutant. The order used here is:

  1. the deposit's sequence matches the enzyme's stored sequence exactly. A variant is a
     different molecule, and an entry named for the wild type is not evidence about it.
  2. failing that, fewest differences from the stored sequence, so a construct differing
     by a purification tag still beats a five-point mutant.
  3. then X-ray over the other methods, and within that the sharpest resolution.

Catalytic knockouts are excluded outright before ranking rather than penalised, for the
same reason they are checked at build time: a S->A construct carries the right name and
the wrong chemistry, and its geometry is a measurement of something that cannot catalyse.
"""

from __future__ import annotations

import json
import urllib.parse
from typing import Dict, List, Optional, Sequence, Tuple

from .. import http
from ..db import connect

UNIPROT_SEARCH = "https://rest.uniprot.org/uniprotkb/search"
RCSB_GRAPHQL = "https://data.rcsb.org/graphql"

# Deposits are fetched in one batched GraphQL call rather than 389 REST calls. The page
# size is bounded because the service returns a 414 on a query string past a few kB.
_BATCH = 40

_QUERY = """query($ids: [String!]!) {
  entries(entry_ids: $ids) {
    rcsb_id
    rcsb_entry_info { resolution_combined experimental_method polymer_entity_count_protein }
    exptl { method }
    polymer_entities {
      entity_poly { pdbx_seq_one_letter_code_can rcsb_entity_polymer_type }
      rcsb_polymer_entity_container_identifiers { auth_asym_ids }
    }
  }
}"""


def pdb_xrefs(accessions: Sequence[str]) -> Dict[str, List[str]]:
    """UniProt accession -> its PDB cross-references, for the accessions that have any."""
    out: Dict[str, List[str]] = {}
    accs = [a.split("-")[0] for a in accessions if a]
    for i in range(0, len(accs), 60):
        chunk = accs[i:i + 60]
        query = " OR ".join(f"accession:{a}" for a in chunk)
        data = http.get_json(
            f"{UNIPROT_SEARCH}?format=json&size=500&fields=accession,xref_pdb"
            f"&query={urllib.parse.quote(query)}")
        for rec in (data or {}).get("results", []):
            ids = [x["id"] for x in rec.get("uniProtKBCrossReferences", [])
                   if x.get("database") == "PDB"]
            if ids:
                out[rec["primaryAccession"]] = ids
    return out


def entry_metadata(pdb_ids: Sequence[str]) -> Dict[str, dict]:
    """Resolution, method and per-chain sequence for each deposit, batched."""
    out: Dict[str, dict] = {}
    ids = sorted({p.upper() for p in pdb_ids})
    for i in range(0, len(ids), _BATCH):
        chunk = ids[i:i + _BATCH]
        resp = http.post_json(RCSB_GRAPHQL, {"query": _QUERY, "variables": {"ids": chunk}})
        if resp.status_code != 200:
            continue
        for e in (resp.json().get("data") or {}).get("entries") or []:
            info = e.get("rcsb_entry_info") or {}
            res = (info.get("resolution_combined") or [None])[0]
            methods = [m.get("method") for m in (e.get("exptl") or []) if m.get("method")]
            seqs = []
            for pe in e.get("polymer_entities") or []:
                poly = pe.get("entity_poly") or {}
                if poly.get("rcsb_entity_polymer_type") != "Protein":
                    continue
                s = (poly.get("pdbx_seq_one_letter_code_can") or "").replace("\n", "")
                if s:
                    seqs.append(s)
            out[e["rcsb_id"]] = {"resolution": res, "methods": methods, "sequences": seqs}
    return out


def _differences(a: str, b: str) -> Tuple[int, int]:
    """(substitutions, residues compared) at the offset that fits the two best.

    Anchoring at position 1 is wrong here and wrong in a way that looks like a result: a
    deposit is usually the mature protein and the stored sequence usually the precursor,
    so a 30-residue signal peptide makes every position disagree and a perfect match
    reports 271 differences. That is the same numbering mismatch `reference.sequence_offset`
    exists to reconcile, one stage earlier -- 7CEF is out by +42 against Cut190 and 4CG1 by
    -40 against TfCut2 -- so the offset is searched for rather than assumed.

    Still ungapped, and that is deliberate: the question is whether this deposit is the
    same molecule as the sequence on record, and a gapped alignment free to insert would
    report an engineered variant as a match. A construct needing an indel scores badly and
    ranks below one that does not, which is the right outcome.

    Coverage is returned alongside so the caller can reject a fragment that agrees
    perfectly over forty residues and says nothing about the rest of the protein.
    """
    best = (len(a) + len(b), 0)
    for off in range(-200, 201):
        lo, hi = max(0, off), min(len(a), len(b) + off)
        if hi - lo < 80:
            continue
        mism = sum(1 for i in range(lo, hi) if a[i] != b[i - off])
        if (mism, -(hi - lo)) < (best[0], -best[1]):
            best = (mism, hi - lo)
    return best


def rank_deposits(sequence: str, ids: Sequence[str], meta: Dict[str, dict],
                  knockouts: Optional[set] = None) -> List[Tuple[str, dict]]:
    """Best deposit first. Entries with no protein chain or no metadata are dropped."""
    scored: List[Tuple[tuple, str, dict]] = []
    for pid in {p.upper() for p in ids}:
        if knockouts and pid in knockouts:
            continue
        m = meta.get(pid)
        if not m or not m["sequences"]:
            continue
        diff, cover = min((_differences(sequence, s) for s in m["sequences"]),
                          key=lambda dc: (dc[0], -dc[1]))
        # A chain covering a third of the protein can match perfectly and still be the
        # wrong evidence, so short overlaps are ranked behind, not treated as exact.
        thin = cover < 0.6 * len(sequence)
        xray = 0 if any("RAY" in (mm or "").upper() for mm in m["methods"]) else 1
        # Missing resolution sorts last within its method rather than first: None would
        # otherwise beat 1.2 A and hand an NMR ensemble the top slot.
        res = m["resolution"] if m["resolution"] is not None else 99.0
        scored.append(((1 if thin else 0, 0 if diff == 0 else 1, diff, xray, res), pid,
                       dict(m, differences=diff, coverage=cover)))
    scored.sort(key=lambda t: t[0])
    return [(pid, m) for _, pid, m in scored]


# How different a deposit may be from the sequence on record and still be treated as that
# enzyme's experimental structure. A handful of positions is an expression construct: a
# cloning scar, a residual tag, the strain the gene was actually taken from. Twenty is a
# different protein, and UniProt cross-references homologues freely -- 8SPK sits under
# P19833 at 20 differences over 282 residues. Past this line the AlphaFold model wins,
# because it was at least computed from THIS enzyme's own sequence, and a crystal structure
# of a near neighbour presented as this enzyme's structure is the worse error: it looks
# like evidence.
MAX_DIFFERENCES = 8
MIN_COVERAGE = 0.60


def link(only: Optional[Sequence[str]] = None, dry_run: bool = False,
         max_differences: int = MAX_DIFFERENCES,
         min_coverage: float = MIN_COVERAGE) -> Dict[str, object]:
    """Populate `pdb_ids_json` for positives that have an accession and no deposit yet.

    Returns a report rather than printing, so a caller can decide what to rebuild: the
    identifiers are only useful once `reference.build(only=...)` reruns for those enzymes
    and the deposit replaces the model on disk.
    """
    with connect() as c:
        sql = ("SELECT enzyme_id, uniprot, sequence FROM characterised_enzymes "
               "WHERE is_positive=1 AND sequence IS NOT NULL "
               "  AND uniprot IS NOT NULL AND uniprot != '' "
               "  AND (pdb_ids_json IS NULL OR pdb_ids_json IN ('', '[]'))")
        params: List[str] = []
        if only:
            sql += f" AND enzyme_id IN ({','.join('?' * len(only))})"
            params = list(only)
        rows = c.execute(sql + " ORDER BY enzyme_id", params).fetchall()

    by_acc: Dict[str, List[tuple]] = {}
    for eid, acc, seq in rows:
        by_acc.setdefault(acc.split("-")[0], []).append((eid, seq))

    xrefs = pdb_xrefs(list(by_acc))
    meta = entry_metadata([p for ids in xrefs.values() for p in ids]) if xrefs else {}

    report: Dict[str, object] = {"considered": len(rows), "with_deposits": 0,
                                 "linked": [], "no_match": [], "too_divergent": []}
    updates: List[Tuple[str, str]] = []
    for acc, entries in by_acc.items():
        ids = xrefs.get(acc)
        if not ids:
            continue
        for eid, seq in entries:
            report["with_deposits"] += 1
            ranked = rank_deposits(seq, ids, meta)
            if not ranked:
                report["no_match"].append(eid)
                continue
            # Only deposits that credibly ARE this enzyme go in, and the list keeps its
            # order: the builder takes the first that returns coordinates, so an entry
            # withdrawn or unreachable falls through to the next best rather than to a
            # prediction.
            ok = [(pid, m) for pid, m in ranked
                  if m["differences"] <= max_differences
                  and m["coverage"] >= min_coverage * len(seq)]
            if not ok:
                worst = ranked[0][1]
                report["too_divergent"].append({
                    "enzyme_id": eid, "uniprot": acc, "pdb_id": ranked[0][0],
                    "differences": worst["differences"], "coverage": worst["coverage"],
                    "seq_length": len(seq)})
                continue
            best_id, best = ok[0]
            report["linked"].append({
                "enzyme_id": eid, "uniprot": acc, "pdb_id": best_id,
                "resolution": best["resolution"], "differences": best["differences"],
                "coverage": best["coverage"], "seq_length": len(seq),
                "alternatives": len(ranked) - 1})
            updates.append((json.dumps([pid for pid, _ in ok[:5]]), eid))

    if updates and not dry_run:
        with connect() as c:
            c.executemany(
                "UPDATE characterised_enzymes SET pdb_ids_json=? WHERE enzyme_id=?", updates)
            c.commit()
    return report
