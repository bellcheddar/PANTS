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
    # UniParc UPI, used ONLY when the UniProtKB accession is dead. UniParc archives
    # sequences and never deletes, so it can still serve a retired entry, but it carries
    # no organism or lineage. Set this and the loader falls back to it, reporting that it
    # did rather than quietly substituting one source for another.
    uniparc: Optional[str] = None


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
        notes="Saccharomonospora viridis AHK190 cutinase, 304 aa. STRAIN AMBIGUITY "
              "RESOLVED 2026-08-05: UniProt carries two 304 aa entries for this organism, "
              "W0TJ64 and C7MVE8 (type strain P101), and length alone cannot separate "
              "them. The Cut190 crystal structures settle it: 4WFI, 4WFJ, 4WFK, 5ZNO, "
              "5ZRQ and 5ZRR all cross-reference W0TJ64, and C7MVE8 has no PDB entry at "
              "all. Structural evidence over a name match.",
        pdb_ids=["4WFI", "4WFJ", "4WFK", "5ZNO", "5ZRQ", "5ZRR"],
    ),
    WildType(
        enzyme_id="BhrPETase", uniprot="A0A2H5Z9R5", uniparc="UPI000CB4D10C",
        family="petase_like", is_positive=True,
        notes="Bacterium HR29 PET hydrolase, 293 aa. Added as a wild type because it is "
              "TurboPETase's parent, and a variant can only be derived from a parent that "
              "lives in WILD_TYPES. It is also a PAZy-measured positive in its own right "
              "(PAZy:17), so it duplicates that row: the same already holds for IsPETase, "
              "LCC and Cut190, and is harmless because evaluation splits by 30% cluster, "
              "which keeps identical sequences on the same side. "
              "NO LIVE UNIPROTKB ENTRY: A0A2H5Z9R5, the accession PAZy still records, is "
              "inactive (DEMERGED to A0ACD6B9U1), and A0ACD6B9U1 is itself DELETED as "
              "'not part of a reference proteome'. The sequence therefore comes from "
              "UniParc, which never deletes, and is confirmed three ways: it is byte-identical "
              "to PAZy:17, it carries an active EMBL WGS cross-reference (GBD22443, "
              "'Poly(Ethylene terephthalate) hydrolase'), and all eight TurboPETase "
              "mutations match it at offset 0.",
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
        notes="ThermoPETase scaffold plus R224Q/N233K. Machine-learning designed. "
              "HAS an experimental structure, contrary to a first look: 7SH6 appears to "
              "carry 19 substitutions against IsPETase, but 14 of those are a heterologous "
              "secretion leader replacing the native signal peptide, and its MATURE chain "
              "(residue 28 on) is identical to FAST-PETase. 8J45 is the same protein again, "
              "expressed in Pichia pastoris, differing only by an EF cloning scar where the "
              "mature sequence begins QT. Comparing whole precursors instead of mature "
              "chains is what hid this.",
        pdb_ids=["7SH6", "8J45"],
    ),
    Variant(
        enzyme_id="FAST-PETase-N212A/K233C/S282C", parent="IsPETase",
        mutations=["S121E", "D186H", "N212A", "R224Q", "N233C", "R280A", "S282C"],
        mutations_confirmed=True,
        reference="Lu et al. 2022, Nature (MutCompute)",
        notes="A disulfide-stabilised FAST-PETase, from PDB 9LMS at 1.71 A. Expressed "
              "against IsPETase rather than against FAST-PETase, because IsPETase is the "
              "root of this lineage and every other member is quoted the same way: it is "
              "FAST-PETase's five with N233K taken on to N233C, plus N212A and S282C. "
              "C233 and C282 form an engineered disulfide across the SAME residue pair "
              "Z1-PETase uses, reached independently by a different group. Verified by "
              "aligning 9LMS's mature chain to FAST-PETase, which differs at exactly the "
              "three positions its deposit title names.",
        pdb_ids=["9LMS"],
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
    Variant(
        enzyme_id="DuraPETase", parent="IsPETase",
        mutations=["S214H", "I168R", "W159H", "S188Q", "R280A",
                   "A180I", "G165A", "Q119Y", "L117F", "T140D"],
        mutations_confirmed=True,
        reference="Cui et al. 2021, ACS Catal.",
        notes="Ten-mutation redesign of IsPETase by the GRAPE strategy, Topt 37 C. All ten "
              "stated parent residues match IsPETase at offset 0, and the count agrees "
              "with the published ten. Ten independent positions agreeing by chance is "
              "~20^-10, so the set is confirmed without needing the supplementary.",
    ),
    Variant(
        enzyme_id="TurboPETase", parent="BhrPETase",
        mutations=["H218S", "F222I", "A209R", "D238K",
                   "A251C", "A281C", "W104L", "F243T"],
        mutations_confirmed=True,
        reference="Zhang et al. 2024, Nat. Commun. 15:1417",
        notes="Eight mutations on BhrPETase, NOT on IsPETase: the parent recorded here "
              "before curation was wrong. Grouped by what they do: substrate-binding "
              "cleft flexibility (H218S/F222I, W104L, F243T), surface charge-charge "
              "optimisation (A209R, D238K) and one engineered disulfide (A251C-A281C). "
              "All eight match BhrPETase at offset 0. doi:10.1038/s41467-024-45662-9",
    ),
    Variant(
        enzyme_id="Z1-PETase", parent="IsPETase",
        mutations=["N37D", "S121E", "R132E", "A171C", "A180V", "P181V", "D186H",
                   "S193C", "R224E", "N233C", "S242T", "N246D", "S282C"],
        mutations_confirmed=True,
        reference="J. Hazard. Mater. 2023",
        notes="Thirteen mutations on IsPETase, including two engineered disulfides "
              "(A171C-S193C, N233C-S282C) and four shared with FAST-PETase "
              "(S121E/D186H/S242T/N246D). Topt 30 C. Confirmed TWICE and independently: "
              "the set applies cleanly to IsPETase at offset 0, and all 13 sites read the "
              "mutant residue in the deposited 8H5K structure, whose sequence differs from "
              "the derived one only by an SHM expression-tag scar and the signal peptide "
              "(zero mismatches in the mature region).",
        pdb_ids=["8H5K"],
    ),
    Variant(
        enzyme_id="DepoPETase", parent="IsPETase",
        mutations=["T88I", "D186H", "D220N", "N233K", "N246D", "R260Y", "S290P"],
        mutations_confirmed=True,
        reference="Cell Rep. Phys. Sci. 2024",
        notes="Seven-mutation IsPETase from flexible-loop directed evolution: Tm +23.3 C "
              "and ~1407-fold more product than wild type. All seven match at offset 0.",
    ),
    Variant(
        enzyme_id="LCC-A2", parent="LCC",
        mutations=["F243I", "D238C", "S283C", "Y127G", "H218Y", "N248D"],
        mutations_confirmed=True,
        reference="Reported relative to LCC-ICCG; 2025 PET-hydrolase review",
        notes="LCC-ICCG plus H218Y/N248D, Topt 78 C. Expressed here against WILD-TYPE LCC "
              "(all six mutations, offset 0) rather than against LCC-ICCG, because a "
              "variant can only be derived from a parent in WILD_TYPES and chaining a "
              "variant onto a variant would hide which residues were actually checked.",
    ),
    # --- recorded, sequence not derived: mutation sets not confirmed during curation ---
    Variant(
        enzyme_id="HotPETase", parent="IsPETase", mutations=[], mutations_confirmed=False,
        reference="Bell et al. 2022, Nat. Catal.",
        notes="Directed-evolution thermostabilised IsPETase, ~21 mutations. Full set not "
              "confirmed here.",
    ),
    Variant(
        enzyme_id="Cut190**SS", parent="Cut190", mutations=[], mutations_confirmed=False,
        reference="Oda/Kawai et al.",
        notes="Disulfide-stabilised Cut190 variant. Full set not confirmed here, and the "
              "parent strain assignment is itself unresolved (see Cut190 notes).",
    ),
]

# PHL7 / PES-H1: DO NOT seed from UniProt. Checked exhaustively 2026-08-05.
#
# PHL7 (Polyester Hydrolase Leipzig 7) is a well-characterised metagenomic polyester
# hydrolase and an obvious candidate parent. It has exactly ONE UniProt entry,
# A0AA82WPD4, and that entry is the **catalysis-deficient S131A mutant** deposited for
# crystallography: the catalytic serine is knocked out. It is 267 aa, the right length,
# named "PHL-7", and differs from the active enzyme at a single position. Seeding from it
# would put a dead enzyme into the positives looking entirely healthy. This is the same
# trap as PDB 7CEH (S176A) rejected during the Cut190**SS work, and length agreement is
# again no defence.
#
# The ACTIVE sequence is already present as PAZy:37, a PAZy-measured positive, so PHL7
# itself is in the training set. It is not promoted to WILD_TYPES because that path
# resolves parents through UniProt and the only accession is the knockout.
#
# Consequence: the PES-H1 L92F/Q94Y variant is not derived here, for want of a loadable
# parent, not for want of a confirmed mutation set. That set was in fact confirmed, and
# resolving it settled a documented literature discrepancy: the same double mutant is
# written L92F/Q94Y in PES-H1 numbering and L93F/Q95Y in PHL7 numbering. find_offset
# returns +1 and 0 respectively against PAZy:37, and both produce an IDENTICAL sequence,
# which is the numbering-offset machinery doing exactly the job it was built for.
#
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


# Variants whose sequence was recovered from a CRYSTAL STRUCTURE rather than derived from
# a mutation list. This is the stronger route where it exists: the PDB SEQRES is the
# construct that was actually expressed, crystallised and assayed, so nothing is applied
# or assumed. The stored sequence is the MATURE construct, not the precursor, which is why
# it is shorter than its parent.
#
# Each was verified by aligning against its parent and checking the substitution count
# against the published one, rather than trusting the name on the PDB entry:
#
#   HotPETase   7QVH  272 aa, 21 substitutions vs IsPETase   (paper reports ~21)
#   Cut190**SS  7CEF  262 aa,  4 substitutions vs Cut190     (S226P/R228S plus two
#                                                             construct-boundary changes)
PDB_DERIVED = {
    "HotPETase": ("7QVH", "IsPETase", 21, "Bell et al. 2022, Nat. Catal."),
    "Cut190**SS": ("7CEF", "Cut190", 4, "Oda/Kawai et al."),
}


# DOIs for the reference strings above, each VERIFIED through Crossref: resolved, and the
# returned type, title, journal and year checked against the citation rather than trusting
# the first search hit.
#
# That check was not ceremony. A bibliographic search returned, for five of these, a DOI
# that resolves to something adjacent but wrong: `10.1021/acscatal.0c05126.s001` is the
# SUPPLEMENTARY FILE for Cui et al., not the paper (the `.sNNN` suffix is the tell);
# `10.26434/chemrxiv-2021-mcjh6` is the ChemRxiv PREPRINT of Bell et al., not the Nature
# Catalysis article; and the top hits for Lu et al. and Tournier et al. were both Faculty
# Opinions RECOMMENDATIONS of the paper, which have their own DOIs and are not the paper.
# All four look right in a citation and send the reader somewhere else.
REFERENCE_DOI = {
    "Cui et al. 2021, ACS Catal.":                 "10.1021/acscatal.0c05126",
    "Son et al. 2019, ACS Catal.":                 "10.1021/acscatal.9b00568",
    "Bell et al. 2022, Nat. Catal.":               "10.1038/s41929-022-00821-3",
    "Lu et al. 2022, Nature (MutCompute)":         "10.1038/s41586-022-04599-z",
    "Tournier et al. 2020, Nature":                "10.1038/s41586-020-2149-4",
    "Austin et al. 2018, PNAS":                    "10.1073/pnas.1718804115",
    "Zhang et al. 2024, Nat. Commun. 15:1417":     "10.1038/s41467-024-45662-9",
    "Cell Rep. Phys. Sci. 2024":                   "10.1016/j.xcrp.2024.102295",
    "J. Hazard. Mater. 2023":                      "10.1016/j.jhazmat.2023.132297",
    "Oda/Kawai et al.":                            "10.1021/acs.biochem.8b00624",
    "Reported relative to LCC-ICCG; 2025 PET-hydrolase review": "10.1002/pro.70282",
    # The wild types' own describing papers, for the pages that have no variant reference.
    "IsPETase":                                    "10.1126/science.aad6359",
}


def count_substitutions(parent: str, construct: str) -> int:
    """Substitutions between a parent and a construct, by GAPPED alignment.

    Must be gapped. A crystallised construct is the mature protein and often carries an
    expression-tag scar, so it is both shorter than the precursor parent and offset from
    it. Counting position-by-position without aligning reported 247 substitutions for
    HotPETase, which really has 21: every residue after the first indel was compared
    against the wrong partner, and the number was wrong in a direction that looks like a
    heavily engineered enzyme rather than like a bug.

    Only columns where BOTH sequences have a residue are counted, so the leading and
    trailing truncations are not miscounted as substitutions.

    biotite is imported lazily: the web venv is deliberately thin and does not have it.
    """
    import biotite.sequence as bseq
    import biotite.sequence.align as balign

    matrix = balign.SubstitutionMatrix.std_protein_matrix()
    alignment = balign.align_optimal(
        bseq.ProteinSequence(parent), bseq.ProteinSequence(construct),
        matrix, gap_penalty=(-10, -1),
    )[0]
    a, b = balign.get_symbols(alignment)
    return sum(1 for x, y in zip(a, b) if x and y and x != y)


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
