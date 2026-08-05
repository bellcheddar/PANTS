"""PAZy: the plastics-active enzyme database.

https://api.pazy.eu/api  ·  Buchholz et al. (2022), Proteins 90(7), 1443-1456.

This is the source the project brief named first, and it is the one that fixes the
labelling problem the rest of the pipeline kept running into.

**Why PAZy entries are evidence where the EC-annotated bulk is not.** The 449 positives
carrying `EC 3.1.1.101` from `ECO:0000256` were assigned by similarity: a sequence model
reproducing them is close to tautological, which is why a head trained on those labels
scored AUC 1.000 and meant nothing. PAZy's inclusion criterion is the opposite. An enzyme
is in it because somebody **measured** activity on a plastic and published it, and every
record carries the DOI.

What that buys, measured rather than hoped for:

| | Evidenced positives | Clusters at 30% | Clusters at 50% |
|---|---|---|---|
| Before | 17 | 5 | 7 |
| With PAZy | 329 | **51** | **71** |

Cluster count is the number that matters, because evaluation splits by cluster and never
by sequence. Five clusters could not support the protocol; fifty-one can.

Note on the `verified` field: every PAZy record has `verified: true`, so it does not
discriminate within the database. It means "curated into PAZy", and PAZy's own criterion
is the experimental one. Do not mistake it for a per-record confidence score.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional

from .. import http
from ..db import connect, now, retry_write
from ..db.manifest import stage_manifest

API = "https://api.pazy.eu/api"
CITATION = "Buchholz et al. 2022, Proteins 90(7):1443-1456"
STAGE = "pazy"

STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")


def fetch_all(page_size: int = 200) -> List[dict]:
    """Every PAZy protein, following pagination."""
    out: List[dict] = []
    url: Optional[str] = f"{API}/proteins/"
    params: Optional[dict] = {"page_size": page_size}
    while url:
        resp = http.get(url, params=params)
        resp.raise_for_status()
        body = resp.json()
        out.extend(body.get("results", []))
        url, params = body.get("next"), None
    return out


def substrates_of(p: dict) -> List[str]:
    return [s.get("abbreviation") or s.get("name") or ""
            for s in (p.get("substrates") or [])]


def _first_doi(p: dict) -> Optional[str]:
    for lit in (p.get("literature") or []):
        doi = lit.get("doi")
        if doi:
            return doi.replace("https://doi.org/", "")
    return None


def load_within_family_negatives(label: str = "v1", identity: float = 0.3
                                 ) -> Dict[str, object]:
    """PAZy enzymes measured on a plastic that is NOT PET, restricted to those sharing a
    sequence cluster with a PET-active enzyme.

    This exists to answer the one question the current evaluation cannot. The head scores
    AUC 0.976 against hard negatives drawn from other alpha/beta-hydrolase families, and
    1.000 against the near misses, but every near miss is a single ESTHER `Cutinase`
    family. Neither result speaks to the question that decides whether PANTS beats a
    homology search: **among enzymes that are already polyesterases, which ones act on
    PET?**

    PAZy can supply that set, because it curates enzymes measured on PA, PUR, PLA, PBAT,
    PHA and others as well as on PET. The restriction to shared clusters is the whole
    point: a nylon amidase or a PLA protease is a different fold doing a different job, and
    including it would be another easy negative inflating the score. Only an enzyme that
    clusters WITH the PET-active set at `identity` is inside the family boundary.

    **The caveat, which must travel with the number.** PAZy records substrates an enzyme
    was SHOWN to degrade. It has no field for a tested-and-inactive result. So "PET not
    listed" means "not reported active on PET", which conflates *inactive* with *never
    assayed*. These are therefore WEAK negatives, and some fraction are false. They are
    stored under their own source_ref so they can always be included or excluded
    deliberately, and any metric computed against them is a lower bound on how well the
    head separates true PET activity, not a clean estimate of it.
    """
    from .. import seqtools

    report: Dict[str, object] = {"added": 0, "candidates": 0, "mixed_clusters": 0,
                                 "skipped": []}

    with stage_manifest(STAGE, label=f"{label}-within-family-negatives",
                        params={"identity": identity}) as m:
        prots = fetch_all()
        clean = [p for p in prots
                 if p.get("amino_acid_sequence")
                 and set(p["amino_acid_sequence"]) <= STANDARD_AA]

        pet = {f"PAZy:{p['id']}" for p in clean if "PET" in substrates_of(p)}
        by_id = {f"PAZy:{p['id']}": p for p in clean}

        import tempfile, pathlib
        tmp = pathlib.Path(tempfile.mkdtemp())
        fasta = seqtools.write_fasta(
            ((k, v["amino_acid_sequence"]) for k, v in by_id.items()), tmp / "pazy.fasta")
        clusters = seqtools.cluster(fasta, min_seq_id=identity)

        grouped: Dict[str, list] = {}
        for member, rep in clusters.items():
            grouped.setdefault(rep, []).append(member)

        keep = []
        for rep, members in grouped.items():
            has_pet = any(mm in pet for mm in members)
            others = [mm for mm in members if mm not in pet]
            if has_pet and others:
                report["mixed_clusters"] = int(report["mixed_clusters"]) + 1
                keep.extend((mm, rep) for mm in others)
        report["candidates"] = len(keep)

        with connect() as conn:
            for eid, rep in keep:
                p = by_id[eid]
                subs = ", ".join(s for s in substrates_of(p) if s)
                doi = _first_doi(p)
                conn.execute(
                    "INSERT INTO characterised_enzymes "
                    "(enzyme_id, uniprot, organism, family, sequence, seq_length, "
                    " is_positive, is_negative, is_near_miss, taxonomy_lineage, "
                    " activity_substrate_notes, source_ref, added_at) "
                    "VALUES (?,?,?,?,?,?,0,0,1,?,?,?,?) "
                    "ON CONFLICT(enzyme_id) DO UPDATE SET "
                    " activity_substrate_notes=excluded.activity_substrate_notes, "
                    " source_ref=excluded.source_ref, is_near_miss=excluded.is_near_miss",
                    (f"{eid}-nonPET", p.get("uniprot_accession") or None,
                     (p.get("organism") or {}).get("scientific_name"), "petase_like",
                     p["amino_acid_sequence"], len(p["amino_acid_sequence"]),
                     (p.get("organism") or {}).get("phylum"),
                     (f"PAZy {p.get('name')}: WITHIN-FAMILY WEAK NEGATIVE. Measured active "
                      f"on {subs}, with PET absent from its measured substrates. Shares a "
                      f"{int(identity*100)}% cluster with PET-active enzymes (cluster "
                      f"representative {rep}), so it is inside the polyesterase family "
                      f"boundary rather than an easy out-of-family negative. "
                      f"CAVEAT: PAZy records only positive substrate associations, so this "
                      f"means 'not reported active on PET', which does not distinguish "
                      f"inactive from never assayed. "
                      + (f"Primary reference doi:{doi}. " if doi else "")
                      + f"Curated set: {CITATION}."),
                     "PAZy-nonPET", now()),
                )
                report["added"] = int(report["added"]) + 1

        m.counts(n_input=len(clean), n_output=int(report["added"]),
                 n_discarded=len(clean) - int(report["added"]))
    return report


def load(substrate: str = "PET", label: str = "v1") -> Dict[str, object]:
    """Write PAZy's enzymes for one substrate into characterised_enzymes.

    Sequences come from PAZy and are stored as-is: non-standard residues are refused
    rather than cleaned, because silently repairing a sequence is how a wrong one enters
    the training set looking healthy.
    """
    report: Dict[str, object] = {"added": 0, "updated": 0, "skipped": [], "substrate": substrate}

    with stage_manifest(STAGE, label=f"{label}-{substrate}") as m:
        prots = fetch_all()
        want = [p for p in prots if substrate in substrates_of(p)]

        with connect() as conn:
            known = {r[0] for r in conn.execute(
                "SELECT uniprot FROM characterised_enzymes WHERE uniprot IS NOT NULL")}

            for p in want:
                seq = p.get("amino_acid_sequence")
                if not seq:
                    report["skipped"].append(f"{p.get('name')}: no sequence")
                    continue
                if not set(seq) <= STANDARD_AA:
                    report["skipped"].append(f"{p.get('name')}: non-standard residues")
                    continue

                acc = p.get("uniprot_accession") or None
                org = (p.get("organism") or {}).get("scientific_name")
                doi = _first_doi(p)
                subs = ", ".join(s for s in substrates_of(p) if s)
                # Keep PAZy's own id in the enzyme_id: it is the stable handle back to the
                # record, and PAZy names are not unique across entries.
                eid = f"PAZy:{p['id']}"
                already = acc in known if acc else False

                conn.execute(
                    "INSERT INTO characterised_enzymes "
                    "(enzyme_id, uniprot, organism, family, sequence, seq_length, "
                    " is_positive, is_negative, is_near_miss, taxonomy_lineage, "
                    " activity_substrate_notes, source_ref, added_at) "
                    "VALUES (?,?,?,?,?,?,1,0,0,?,?,?,?) "
                    "ON CONFLICT(enzyme_id) DO UPDATE SET sequence=excluded.sequence, "
                    " seq_length=excluded.seq_length, "
                    " activity_substrate_notes=excluded.activity_substrate_notes",
                    (eid, acc, org, "petase_like", seq, len(seq),
                     (p.get("organism") or {}).get("phylum"),
                     (f"PAZy {p.get('name')}: measured activity on {subs}. "
                      f"In PAZy because activity was MEASURED and published, not because "
                      f"of sequence similarity. "
                      + (f"Primary reference doi:{doi}. " if doi else "")
                      + f"Curated set: {CITATION}."),
                     "PAZy-measured", now()),
                )
                report["updated" if already else "added"] += 1

        def _src() -> None:
            with connect() as conn:
                conn.execute(
                    "INSERT INTO data_sources(name, version, retrieved_at, n_records, "
                    "license, source_url) VALUES (?,?,?,?,?,?) "
                    "ON CONFLICT(name) DO UPDATE SET retrieved_at=excluded.retrieved_at, "
                    "n_records=excluded.n_records",
                    ("PAZy", CITATION, now(), len(want), "see publication", API))
        retry_write(_src)

        m.counts(n_input=len(prots), n_output=int(report["added"]) + int(report["updated"]),
                 n_discarded=len(report["skipped"]))
    return report
