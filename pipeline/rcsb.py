"""Minimal RCSB PDB REST access: the deposited sequence for a structure.

Used by the seed loader for variants whose sequence is recovered from a crystal structure
rather than derived from a mutation list. The PDB SEQRES is the construct that was
actually expressed, crystallised and assayed, which makes it stronger evidence than a
mutation list transcribed out of a paper, provided the entry really is the enzyme claimed.

That proviso is not theoretical. Two entries were rejected during curation for looking
right and being wrong: 7CEH is a Cut190 variant carrying S176A, the catalytic serine
knocked out, and the sole UniProt entry for PHL7 is likewise an S131A knockout. Both are
the correct length and carry the expected name. So the caller must verify the substitution
count against the published one rather than trusting the title, which is why this module
only fetches and never decides.
"""

from __future__ import annotations

from typing import Optional

from . import http

ENTITY_URL = "https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity}"


def entity_sequence(pdb_id: str, entity: int = 1) -> Optional[str]:
    """One-letter canonical sequence of a polymer entity, or None if unavailable.

    The `_can` form is requested so modified residues come back as their standard parent
    rather than as X, which would otherwise trip the non-standard-residue guard on an
    otherwise perfectly good construct.
    """
    obj = http.get_json(ENTITY_URL.format(pdb_id=pdb_id.upper(), entity=entity))
    if not obj:
        return None
    seq = (obj.get("entity_poly") or {}).get("pdbx_seq_one_letter_code_can")
    return seq.replace("\n", "").strip() if seq else None
