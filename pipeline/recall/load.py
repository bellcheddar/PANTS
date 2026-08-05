"""Populate `characterised_enzymes` with the wild-type seeds and derived variants.

Idempotent: re-running replaces rows by enzyme_id rather than duplicating them, so the
curation can be re-run whenever a mutation set is confirmed.
"""

from __future__ import annotations

import json
from typing import Dict

from .. import config, rcsb, uniprot
from ..db import connect, now, retry_write
from ..db.manifest import stage_manifest
from . import seeds

STAGE = "seeds"


def _upsert(conn, **row) -> None:
    cols = ", ".join(row)
    marks = ", ".join("?" * len(row))
    updates = ", ".join(f"{c}=excluded.{c}" for c in row if c != "enzyme_id")
    conn.execute(
        f"INSERT INTO characterised_enzymes ({cols}) VALUES ({marks}) "
        f"ON CONFLICT(enzyme_id) DO UPDATE SET {updates}",
        tuple(row.values()),
    )


def load_family_positives(label: str = "v1", max_scan: int = 40000,
                          length_min: int = 200, length_max: int = 500) -> Dict[str, object]:
    """Add the whole ESTHER Polyesterase-lipase-cutinase family as annotation-only positives.

    Why this exists: the nine hand-curated positives (IsPETase and its variants, LCC and
    LCC-ICCG, TfCut2, Cut190) collapse into ONE cluster at both 30% and 50% identity. By
    spec section 8's own rule (split by cluster, never by sequence) that is a single
    independent example, so no cluster-split evaluation is possible and any
    cross-validation over them is pure leakage.

    The family harvest brings the positive set to something with real cluster diversity.
    The trade is label quality: these carry ESTHER family annotation, not measured PET
    activity, so they are marked source_ref='ESTHER-family' and must be reported
    separately from the characterised subset (spec section 8's last bullet).
    """
    from ..negatives import esther

    report: Dict[str, object] = {"added": 0, "skipped": 0}
    with stage_manifest("positives_family", label=label,
                        params={"length_min": length_min, "length_max": length_max}) as m:
        hits = [
            h for h in esther.stream(
                f"database:esther AND length:[{length_min} TO {length_max}]",
                max_results=max_scan)
            if h.family in esther.POSITIVE_FAMILIES and h.is_clean
        ]
        with connect() as conn:
            known = {r[0] for r in conn.execute(
                "SELECT uniprot FROM characterised_enzymes WHERE uniprot IS NOT NULL")}
            for h in hits:
                if h.accession in known:
                    report["skipped"] += 1      # already present as a curated positive
                    continue
                _upsert(
                    conn,
                    enzyme_id=f"PLC:{h.accession}", uniprot=h.accession,
                    organism=h.organism, family="petase_like",
                    sequence=h.sequence, seq_length=h.length,
                    is_positive=1, is_negative=0, is_near_miss=0,
                    esther_family=h.family, taxonomy_lineage=h.lineage,
                    activity_substrate_notes=(
                        f"ANNOTATION ONLY: ESTHER family Polyesterase-lipase-cutinase. "
                        f"No measured PET activity. UniProt evidence: "
                        f"{h.protein_existence or 'unknown'}"
                        f"{', Swiss-Prot reviewed' if h.reviewed else ''}. "
                        f"Report separately from the characterised subset."),
                    # Evidence level is the axis spec section 8's last bullet needs:
                    # 'Evidence at protein level' is a materially stronger label than
                    # 'Predicted', and the two must not be pooled in a reported metric.
                    source_ref=("ESTHER-family-protein-evidence"
                                if (h.protein_existence or "").startswith("1")
                                else "ESTHER-family-predicted"),
                    added_at=now(),
                )
                report["added"] += 1
        m.counts(n_input=len(hits), n_output=int(report["added"]),
                 n_discarded=int(report["skipped"]))
    return report


def load_reference_set(label: str = "v1") -> Dict[str, object]:
    """Fetch every wild type from UniProt, derive the confirmed variants, write both.

    Returns a report dict; nothing is printed here so the CLI owns presentation.
    """
    report: Dict[str, object] = {"wild_types": [], "variants": [], "problems": []}

    with stage_manifest(STAGE, label=label) as m:
        wt_seqs: Dict[str, str] = {}

        with connect() as conn:
            # --- wild types: real UniProt entries ---
            for w in seeds.WILD_TYPES:
                entry = uniprot.fetch(w.uniprot)
                # A retired accession comes back as an entry with an EMPTY sequence, not
                # as a 404, so "not found" is not the only failure to handle.
                if (entry is None or not entry.sequence) and w.uniparc:
                    entry = uniprot.fetch_uniparc(w.uniparc)
                    report["problems"].append(
                        f"{w.enzyme_id}: UniProtKB {w.uniprot} is dead, fell back to "
                        f"UniParc {w.uniparc} (no organism or lineage from that source)"
                    )
                if entry is None:
                    report["problems"].append(f"{w.enzyme_id}: UniProt {w.uniprot} not found")
                    continue
                if not entry.is_plausible_protein:
                    report["problems"].append(
                        f"{w.enzyme_id}: sequence carries non-standard residues, refused"
                    )
                    continue

                # A recorded triad must actually read Ser/Asp/His in the fetched sequence.
                # If it does not, the accession or the numbering is wrong, and everything
                # downstream (profiles, geometry) would inherit the error.
                if w.triad:
                    got = "".join(entry.sequence[p - 1] for p in w.triad)
                    if sorted(got) != sorted("SDH"):
                        report["problems"].append(
                            f"{w.enzyme_id}: triad {w.triad} reads {got}, not Ser/Asp/His"
                        )
                        continue

                wt_seqs[w.enzyme_id] = entry.sequence
                _upsert(
                    conn,
                    enzyme_id=w.enzyme_id, uniprot=w.uniprot,
                    pdb_ids_json=json.dumps(w.pdb_ids), organism=entry.organism,
                    family=w.family, sequence=entry.sequence, seq_length=entry.length,
                    is_positive=int(w.is_positive), is_negative=0, is_near_miss=0,
                    taxonomy_lineage=entry.lineage,
                    activity_substrate_notes=w.notes, source_ref="UniProt", added_at=now(),
                )
                report["wild_types"].append((w.enzyme_id, entry.length, entry.organism))

            # --- variants: derived from a parent, never pasted ---
            for v in seeds.VARIANTS:
                parent_seq = wt_seqs.get(v.parent)
                if parent_seq is None:
                    report["problems"].append(f"{v.enzyme_id}: parent {v.parent} unavailable")
                    continue

                seq, offset, status = seeds.derive_variant(parent_seq, v)
                if status.startswith("failed"):
                    report["problems"].append(f"{v.enzyme_id}: {status}")

                notes = v.notes
                if status == "unconfirmed":
                    notes = ("MUTATION SET NOT CONFIRMED: no sequence derived, excluded "
                             "from training. " + notes)
                elif status == "derived" and offset:
                    notes = f"Derived at numbering offset {offset:+d}. " + notes

                _upsert(
                    conn,
                    enzyme_id=v.enzyme_id, uniprot=None,
                    pdb_ids_json=json.dumps(v.pdb_ids),
                    organism=None,
                    family=next((w.family for w in seeds.WILD_TYPES
                                 if w.enzyme_id == v.parent), None),
                    sequence=seq, seq_length=len(seq) if seq else None,
                    # Only a sequence-resolved variant counts as a training positive.
                    is_positive=int(seq is not None), is_negative=0, is_near_miss=0,
                    matched_positive_id=v.parent,
                    activity_substrate_notes=notes,
                    source_ref=v.reference, added_at=now(),
                )
                report["variants"].append((v.enzyme_id, status, offset,
                                           len(seq) if seq else None))

            # --- variants recovered from a crystal structure rather than a mutation list ---
            # This runs INSIDE the same load, not as a follow-up script. It was a one-off
            # script once, and re-running the loader silently reset HotPETase and
            # Cut190**SS to sequence-less rows, because the VARIANTS entries that write
            # them first are deliberately unconfirmed. A stage that must be remembered
            # separately is a stage that will eventually be forgotten.
            for eid, (pdb_id, parent, expected, ref) in seeds.PDB_DERIVED.items():
                parent_seq = wt_seqs.get(parent)
                seq = rcsb.entity_sequence(pdb_id)
                if not seq or parent_seq is None:
                    report["problems"].append(f"{eid}: {pdb_id} sequence unavailable")
                    continue
                n_sub = seeds.count_substitutions(parent_seq, seq)
                # The published substitution count is the check. A wrong PDB entry (a
                # catalytic knockout, say) shows up here as a count that does not agree.
                if n_sub != expected:
                    report["problems"].append(
                        f"{eid}: {pdb_id} gives {n_sub} substitutions vs {parent}, "
                        f"expected {expected}; NOT stored"
                    )
                    continue
                _upsert(
                    conn, enzyme_id=eid, uniprot=None, pdb_ids_json=json.dumps([pdb_id]),
                    organism=None,
                    family=next((w.family for w in seeds.WILD_TYPES
                                 if w.enzyme_id == parent), None),
                    sequence=seq, seq_length=len(seq),
                    is_positive=1, is_negative=0, is_near_miss=0,
                    matched_positive_id=parent,
                    activity_substrate_notes=(
                        f"Sequence taken from the deposited construct {pdb_id}, not derived "
                        f"from a mutation list: the PDB SEQRES is what was actually "
                        f"expressed, crystallised and assayed. {n_sub} substitutions vs "
                        f"{parent}, matching the published {expected}. Mature construct, "
                        f"so shorter than the precursor parent."),
                    source_ref="PDB-construct", added_at=now(),
                )
                report["variants"].append((eid, f"pdb:{pdb_id}", None, len(seq)))

        n_derived = sum(1 for _, s, _, _ in report["variants"]
                        if s == "derived" or s.startswith("pdb:"))
        m.counts(n_input=len(seeds.WILD_TYPES) + len(seeds.VARIANTS),
                 n_output=len(report["wild_types"]) + n_derived,
                 n_discarded=sum(1 for _, s, _, _ in report["variants"] if s != "derived"))

        def _src() -> None:
            with connect() as conn:
                conn.execute(
                    "INSERT INTO data_sources(name, version, retrieved_at, n_records, "
                    "license, source_url) VALUES (?,?,?,?,?,?) "
                    "ON CONFLICT(name) DO UPDATE SET retrieved_at=excluded.retrieved_at, "
                    "n_records=excluded.n_records",
                    ("UniProt", "REST", now(), len(wt_seqs),
                     "CC BY 4.0", config.UNIPROT_REST_URL),
                )
        retry_write(_src)

    return report
