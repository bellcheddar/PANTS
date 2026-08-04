"""Catalytic-site anchors, taken from UniProt rather than assumed.

Each profile HMM needs a reference sequence whose catalytic triad positions are known, so
that aligning a candidate to that profile and reading the reference's columns means
something. Hardcoding one reference (IsPETase) only works for profiles that IsPETase
aligns to, which is why the classic `Cutinase` family scored 0/111 on the first pass: the
near misses never aligned to a Polyesterase-lipase-cutinase profile at all.

UniProt annotates `Active site` features with `Nucleophile` and `Charge relay system`
descriptions, so the anchor can be looked up per family instead of guessed.

This was cross-checked against the alignment method before being adopted. Aligning to a
pooled profile and reading IsPETase's verified S160/D206/H237 columns predicted
LCC S165/D210/H242 and TfCut2 S170/D216/H248; UniProt's independent curated annotations
give exactly those numbers. Two methods, three enzymes, no disagreement.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple

from . import seeds
from .. import http, uniprot

ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"

NUCLEOPHILE_HINTS = ("nucleophile", "acyl-ester intermediate")
RELAY_HINTS = ("charge relay",)


@dataclass
class Anchor:
    """A reference sequence with known catalytic positions, in its own 1-based numbering."""
    sequence_id: str
    accession: Optional[str]
    sequence: str
    ser: int
    asp: int
    his: int
    source: str          # 'uniprot' | 'propagated' | 'manual'

    def residues(self) -> str:
        return "".join(self.sequence[p - 1] for p in (self.ser, self.asp, self.his))

    @property
    def is_consistent(self) -> bool:
        """The anchor must actually read Ser/Asp/His in its own sequence."""
        return sorted(self.residues()) == sorted("SDH")


def fetch_active_sites(accession: str) -> Tuple[Optional[str], List[Tuple[int, str]]]:
    """(sequence, [(position, description)]) for a UniProt accession."""
    obj = http.get_json(ENTRY_URL.format(accession=accession))
    if not obj:
        return None, []
    seq = (obj.get("sequence") or {}).get("value", "")
    sites: List[Tuple[int, str]] = []
    for f in obj.get("features", []):
        if f.get("type") != "Active site":
            continue
        pos = (f.get("location", {}).get("start", {}) or {}).get("value")
        if pos:
            sites.append((int(pos), (f.get("description") or "").lower()))
    return seq, sites


def anchor_from_uniprot(sequence_id: str, accession: str) -> Optional[Anchor]:
    """Build an anchor from UniProt's own active-site annotation, or None if unusable.

    The triad is assigned by residue identity rather than by trusting the description
    strings: descriptions vary ('Nucleophile', 'Acyl-ester intermediate') and the charge
    relay pair does not say which member is the Asp and which the His.
    """
    seq, sites = fetch_active_sites(accession)
    if not seq or len(sites) < 3:
        return None

    ser = asp = his = None
    for pos, _desc in sites:
        if not 1 <= pos <= len(seq):
            return None
        aa = seq[pos - 1]
        if aa == "S" and ser is None:
            ser = pos
        elif aa == "D" and asp is None:
            asp = pos
        elif aa == "H" and his is None:
            his = pos

    if None in (ser, asp, his):
        return None
    a = Anchor(sequence_id=sequence_id, accession=accession, sequence=seq,
               ser=ser, asp=asp, his=his, source="uniprot")
    return a if a.is_consistent else None


def find_anchor_for_cluster(members: List[str], sequences: Dict[str, str],
                            accessions: Dict[str, Optional[str]]) -> Optional[Anchor]:
    """First member of a cluster that carries a usable UniProt active-site annotation.

    Members are tried in the order given, so callers should pass the most trustworthy
    first (curated wild types before bulk family entries).
    """
    for m in members:
        acc = accessions.get(m)
        if not acc:
            continue
        a = anchor_from_uniprot(m, acc)
        if a is not None:
            # Guard against UniProt's sequence differing from the one we clustered on.
            if sequences.get(m) and a.sequence != sequences[m]:
                continue
            return a
    return None


# IsPETase's triad, verified twice: read directly off the fetched sequence during Phase 1
# curation, and independently annotated as Active site by UniProt. Used as the fallback
# anchor and as the reference in tests.
ISPETASE = ("IsPETase", "A0A0K8P6T7", 160, 206, 237)
