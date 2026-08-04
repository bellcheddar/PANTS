"""The characterised seed set: wild-type polyesterases and the engineered variants.

Two kinds of entry, handled differently on purpose:

  WILD_TYPES are real UniProt entries. Their sequences are fetched, never typed in.

  VARIANTS (FAST-PETase, LCC-ICCG, DuraPETase and the rest) are engineered mutants. They
  have no UniProt accession of their own, so each is stored as a parent plus a mutation
  list, and the sequence is derived by applying the mutations to the parent.

Deriving rather than pasting is what makes the mutation sets checkable. `apply_mutations`
refuses to substitute unless the parent already carries the residue the mutation names:
a wrong position, or literature numbering that counts from the mature protein while the
parent is a precursor, fails loudly here instead of silently producing a plausible but
wrong sequence that would then poison the HMM profiles and every embedding downstream.

Where the complete mutation set could not be confirmed, the variant is recorded with
mutations_confirmed=False and NO sequence. A partial mutation set yields a wrong sequence,
which is worse than an honest gap: such entries appear in the catalogue as known but not
yet sequence-resolved, and are excluded from training.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

_MUT = re.compile(r"^([ACDEFGHIKLMNPQRSTVWY])(\d+)([ACDEFGHIKLMNPQRSTVWY])$")


class MutationError(ValueError):
    """A mutation does not match the parent sequence at the position it names."""


@dataclass
class WildType:
    enzyme_id: str
    uniprot: str
    family: str            # petase_like | mhetase_like | cutinase | lipase | carboxylesterase
    is_positive: bool
    notes: str = ""
    pdb_ids: List[str] = field(default_factory=list)
    # 1-based positions of the catalytic triad in the UNIPROT PRECURSOR numbering,
    # where known and verified. Used to sanity check the recall stage's triad detection.
    triad: Optional[Tuple[int, int, int]] = None


@dataclass
class Variant:
    enzyme_id: str
    parent: str            # enzyme_id of the WildType it derives from
    mutations: List[str]
    mutations_confirmed: bool
    reference: str
    notes: str = ""
    pdb_ids: List[str] = field(default_factory=list)


# --------------------------------------------------------------------------------------
# Wild types. Accessions verified against UniProt: organism, length and protein name all
# checked, because an accession that merely resolves is not an accession that is correct.
# Q6A0I4 was initially mislabelled Cut190 during this curation; it is TfCut2.
# --------------------------------------------------------------------------------------
WILD_TYPES: List[WildType] = [
    WildType(
        enzyme_id="IsPETase", uniprot="A0A0K8P6T7", family="petase_like", is_positive=True,
        # UniProt now lists the organism as Piscinibacter sakaiensis (reclassified from
        # Ideonella sakaiensis, under which name the entire literature is written).
        notes="Ideonella/Piscinibacter sakaiensis 201-F6. The reference PETase. 290 aa "
              "precursor; literature mutation numbering matches this precursor directly "
              "(verified: S160/D206/H237 triad, mobile W185).",
        pdb_ids=["6EQE", "5XJH", "6ANE"], triad=(160, 206, 237),
    ),
    WildType(
        enzyme_id="MHETase", uniprot="A0A0K8P8E7", family="mhetase_like", is_positive=True,
        notes="Ideonella/Piscinibacter sakaiensis. Tannase/feruloyl esterase group with a "
              "lid domain, NOT the PETase family: a PETase-seeded profile search never "
              "reaches it (spec section 2, point 5). 603 aa.",
        pdb_ids=["6QGA", "6QZ3"],
    ),
    WildType(
        enzyme_id="LCC", uniprot="G9BY57", family="cutinase", is_positive=True,
        notes="Leaf-branch compost cutinase, metagenome-derived. Parent of the LCC-ICCG "
              "industrial variant. 293 aa.",
        pdb_ids=["4EB0"],
    ),
    WildType(
        enzyme_id="TfCut2", uniprot="Q6A0I4", family="cutinase", is_positive=True,
        notes="Thermobifida fusca cutinase 2. Thermophilic, industrial lineage. 301 aa.",
        pdb_ids=["4CG1"],
    ),
    WildType(
        enzyme_id="Cut190", uniprot="W0TJ64", family="cutinase", is_positive=True,
        notes="Saccharomonospora viridis cutinase, 304 aa. STRAIN CAVEAT: Cut190 is from "
              "S. viridis AHK190; UniProt also carries C7MVE8 (type strain P101) at the "
              "same length. Confirm the strain against the Cut190 papers before this seed "
              "is used for anything load-bearing.",
        pdb_ids=["4WFI", "4WFJ"],
    ),
]

# --------------------------------------------------------------------------------------
# Engineered variants.
#
# mutations_confirmed=True means the COMPLETE mutation set is believed correct and is
# checked against the parent by apply_mutations(). False means the variant is real and
# well known but its full mutation set was not confirmed during curation: it is recorded
# for the catalogue, carries no derived sequence, and is excluded from training.
# --------------------------------------------------------------------------------------
VARIANTS: List[Variant] = [
    Variant(
        enzyme_id="ThermoPETase", parent="IsPETase",
        mutations=["S121E", "D186H", "R280A"], mutations_confirmed=True,
        reference="Son et al. 2019, ACS Catal.",
        notes="Thermostabilised IsPETase, the scaffold FAST-PETase was built on.",
    ),
    Variant(
        enzyme_id="FAST-PETase", parent="IsPETase",
        mutations=["S121E", "D186H", "R224Q", "N233K", "R280A"], mutations_confirmed=True,
        reference="Lu et al. 2022, Nature (MutCompute)",
        notes="ThermoPETase scaffold plus R224Q/N233K. Machine-learning designed.",
    ),
    Variant(
        enzyme_id="IsPETase-W159H/S238F", parent="IsPETase",
        mutations=["W159H", "S238F"], mutations_confirmed=True,
        reference="Austin et al. 2018, PNAS",
        notes="The narrowed-cleft double mutant: the experiment that showed the PETase "
              "cleft is wider than a cutinase's and that narrowing it changes activity.",
        pdb_ids=["6EQF"],
    ),
    Variant(
        enzyme_id="LCC-ICCG", parent="LCC",
        mutations=["F243I", "D238C", "S283C", "Y127G"], mutations_confirmed=True,
        reference="Tournier et al. 2020, Nature",
        notes="The industrial benchmark: D238C/S283C disulfide plus F243I/Y127G. Numbering "
              "is checked against the parent by apply_mutations; if the literature counts "
              "from the mature protein the offset is reported rather than guessed at.",
        pdb_ids=["6THT"],
    ),
    # --- recorded, sequence not derived: mutation sets not confirmed during curation ---
    Variant(
        enzyme_id="DuraPETase", parent="IsPETase", mutations=[], mutations_confirmed=False,
        reference="Cui et al. 2021, ACS Catal.",
        notes="Ten-mutation redesign of IsPETase (GRAPE strategy). Full set not confirmed "
              "here; needs the paper's supplementary table before a sequence is derived.",
    ),
    Variant(
        enzyme_id="HotPETase", parent="IsPETase", mutations=[], mutations_confirmed=False,
        reference="Bell et al. 2022, Nat. Catal.",
        notes="Directed-evolution thermostabilised IsPETase, ~21 mutations. Full set not "
              "confirmed here.",
    ),
    Variant(
        enzyme_id="TurboPETase", parent="IsPETase", mutations=[], mutations_confirmed=False,
        reference="Zhang et al. 2024",
        notes="Full mutation set not confirmed here.",
    ),
    Variant(
        enzyme_id="Z1-PETase", parent="IsPETase", mutations=[], mutations_confirmed=False,
        reference="Literature, to confirm",
        notes="Full mutation set not confirmed here.",
    ),
    Variant(
        enzyme_id="Cut190**SS", parent="Cut190", mutations=[], mutations_confirmed=False,
        reference="Oda/Kawai et al.",
        notes="Disulfide-stabilised Cut190 variant. Full set not confirmed here, and the "
              "parent strain assignment is itself unresolved (see Cut190 notes).",
    ),
]

# HGMP01 is named in spec section 1 but is metagenome-derived rather than a variant of a
# characterised parent, so it does not belong in VARIANTS.
#
# Checked exhaustively 2026-08-04. HGMP01's sequence is NOT obtainable by any programmatic
# route. Every avenue tried and its result:
#
#   UniProt, 4 query forms (name / gene / protein_name / gut+EC)   0 hits
#   NCBI protein                                                   0 hits
#   NCBI nuccore                                                   0 hits
#   NCBI elink from PMID 39551294                                  no sequence records,
#                                                                  only pubmed cross-links
#   PubMed record                                                  no accessions, no DAS
#   ScienceDirect full text                                        HTTP 403, subscription
#   Europe PMC                                                     isOpenAccess N, inEPMC N,
#                                                                  hasSuppl N, no data
#                                                                  availability statement
#   MGnify protein API                                             endpoint does not exist
#                                                                  (404; the API exposes
#                                                                  studies/assemblies only)
#
# The paper is Int J Biol Macromol 2024, PMID 39551294, doi 10.1016/j.ijbiomac.2024.137732.
# The sequence is in a supplementary PDF behind a subscription. Retrieving it needs
# institutional access, which is a human step, not an automatable one.
#
# It is deliberately NOT reconstructed: a fabricated seed would poison the profiles and
# every score downstream, and would be undetectable afterwards.
#
# That turns out to be an advantage. The same paper reports HGMP01-LIKE genes as widely
# distributed across the human gut microbiome, so running recall over gut assemblies
# WITHOUT the seed is a blind test of whether the pipeline finds them unaided. Seeding
# with the sequence would have made that test circular.
#
# TODO: pull the sequence from the paper's supplementary by hand for a direct comparison
# against whatever the gut recall returns.


def parse_mutation(mut: str) -> Tuple[str, int, str]:
    m = _MUT.match(mut.strip())
    if not m:
        raise MutationError(f"unparseable mutation {mut!r}: expected e.g. 'S121E'")
    return m.group(1), int(m.group(2)), m.group(3)


def find_offset(sequence: str, mutations: List[str], search: range = range(-60, 61)) -> Optional[int]:
    """Smallest offset making every mutation's stated parent residue match.

    Literature numbering often counts from the mature protein while UniProt stores the
    precursor, so the whole mutation set is shifted by the signal peptide length. Rather
    than guessing that length, this determines it: an offset that satisfies every mutation
    at once is strong evidence, whereas a single mutation matching by chance is not.
    """
    for off in sorted(search, key=abs):
        if all(
            0 <= (pos + off - 1) < len(sequence) and sequence[pos + off - 1] == wt
            for wt, pos, _ in (parse_mutation(m) for m in mutations)
        ):
            return off
    return None


def apply_mutations(sequence: str, mutations: List[str], offset: int = 0) -> str:
    """Apply substitutions, refusing any whose stated parent residue does not match."""
    seq = list(sequence)
    for mut in mutations:
        wt, pos, new = parse_mutation(mut)
        idx = pos + offset - 1
        if not 0 <= idx < len(seq):
            raise MutationError(f"{mut}: position {pos} (offset {offset}) outside 1..{len(seq)}")
        if seq[idx] != wt:
            raise MutationError(
                f"{mut}: parent has {seq[idx]} at position {pos} (offset {offset}), not {wt}"
            )
        seq[idx] = new
    return "".join(seq)


def derive_variant(parent_sequence: str, variant: Variant) -> Tuple[Optional[str], Optional[int], str]:
    """Return (sequence, offset_used, status) for a variant.

    status is one of: 'derived', 'unconfirmed', 'failed:<reason>'. A failure never raises:
    the curation run reports every problem at once rather than dying on the first.
    """
    if not variant.mutations_confirmed or not variant.mutations:
        return None, None, "unconfirmed"
    try:
        off = find_offset(parent_sequence, variant.mutations)
        if off is None:
            return None, None, "failed:no offset satisfies all mutations"
        return apply_mutations(parent_sequence, variant.mutations, off), off, "derived"
    except MutationError as exc:
        return None, None, f"failed:{exc}"
