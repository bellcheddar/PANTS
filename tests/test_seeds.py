"""Tests for the variant-derivation machinery.

These matter more than they look. A wrong mutation set produces a sequence that is still a
valid protein, still folds in ESMFold, still embeds in ESM-2 and still trains: the error
would never surface as a crash, only as quietly degraded scores. The residue-match check
is the only thing standing between a typo and a poisoned seed set.
"""

from __future__ import annotations

import pytest

from pipeline.recall import seeds
from pipeline.recall.seeds import MutationError, Variant


# A toy parent: position 1..10 = A C D E F G H I K L
PARENT = "ACDEFGHIKL"


def test_parse_mutation_roundtrip():
    assert seeds.parse_mutation("S121E") == ("S", 121, "E")


@pytest.mark.parametrize("bad", ["", "S121", "121E", "X121E", "S121X", "SE121", "s121e "])
def test_parse_mutation_rejects_malformed(bad):
    with pytest.raises(MutationError):
        seeds.parse_mutation(bad)


def test_apply_mutations_substitutes_correctly():
    assert seeds.apply_mutations(PARENT, ["A1M"]) == "MCDEFGHIKL"
    assert seeds.apply_mutations(PARENT, ["A1M", "L10V"]) == "MCDEFGHIKV"


def test_apply_mutations_refuses_wrong_parent_residue():
    """The whole point: a mutation naming a residue the parent does not have is refused."""
    with pytest.raises(MutationError, match=r"parent has A at position 1 \(offset 0\), not S"):
        seeds.apply_mutations(PARENT, ["S1E"])


def test_apply_mutations_refuses_out_of_range():
    with pytest.raises(MutationError, match="outside"):
        seeds.apply_mutations(PARENT, ["A99M"])


def test_apply_mutations_leaves_parent_unmodified():
    seeds.apply_mutations(PARENT, ["A1M"])
    assert PARENT == "ACDEFGHIKL"


def test_find_offset_detects_mature_protein_numbering():
    """Literature often numbers from the mature protein while UniProt stores the precursor.
    A consistent shift across every mutation is what identifies the offset."""
    # PARENT is A C D E F G H I K L, so E is at position 4 and G at position 6. Written in
    # numbering shifted by -3, those become E1 and G3, and only offset +3 satisfies both.
    assert seeds.find_offset(PARENT, ["E1F", "G3H"]) == 3


def test_find_offset_returns_none_when_no_shift_works():
    assert seeds.find_offset(PARENT, ["W1A", "W2A"]) is None


def test_find_offset_prefers_zero_when_zero_works():
    assert seeds.find_offset(PARENT, ["A1M"]) == 0


def test_derive_variant_unconfirmed_yields_no_sequence():
    """An unconfirmed mutation set must produce NO sequence: a partial set would give a
    wrong one, which is worse than an honest gap."""
    v = Variant(enzyme_id="X", parent="P", mutations=[], mutations_confirmed=False,
                reference="")
    seq, off, status = seeds.derive_variant(PARENT, v)
    assert (seq, off, status) == (None, None, "unconfirmed")


def test_derive_variant_reports_failure_without_raising():
    v = Variant(enzyme_id="X", parent="P", mutations=["W1A"], mutations_confirmed=True,
                reference="")
    seq, _, status = seeds.derive_variant(PARENT, v)
    assert seq is None and status.startswith("failed")


def test_every_confirmed_variant_has_mutations():
    """mutations_confirmed=True with an empty list would silently mark a variant as
    sequence-resolved while deriving nothing."""
    for v in seeds.VARIANTS:
        if v.mutations_confirmed:
            assert v.mutations, f"{v.enzyme_id} is confirmed but lists no mutations"


def test_every_variant_parent_exists():
    known = {w.enzyme_id for w in seeds.WILD_TYPES}
    for v in seeds.VARIANTS:
        assert v.parent in known, f"{v.enzyme_id} names unknown parent {v.parent}"


def test_wild_type_ids_and_accessions_are_unique():
    ids = [w.enzyme_id for w in seeds.WILD_TYPES] + [v.enzyme_id for v in seeds.VARIANTS]
    assert len(ids) == len(set(ids))
    accs = [w.uniprot for w in seeds.WILD_TYPES]
    assert len(accs) == len(set(accs))


# --------------------------------------------------------------------------------------
# count_substitutions: the gapped-alignment requirement
# --------------------------------------------------------------------------------------

def test_count_substitutions_ignores_terminal_truncation():
    """A mature construct is shorter at both ends than its precursor parent. Those missing
    residues are not substitutions, and counting them as such is what made an ungapped
    comparison report 247 differences for HotPETase, which has 21."""
    parent = "MKKLLAAWQTPYNARGPDPTAASLEASAG"
    construct = parent[6:-3]                       # trimmed both ends, no substitutions
    assert seeds.count_substitutions(parent, construct) == 0


def test_count_substitutions_counts_real_substitutions():
    parent = "MKKLLAAWQTPYNARGPDPTAASLEASAG"
    construct = parent[:10] + "W" + parent[11:]    # exactly one change
    assert construct != parent
    assert seeds.count_substitutions(parent, construct) == 1


def test_count_substitutions_survives_an_internal_deletion():
    """The failure mode being guarded: without gapping, every residue after an indel is
    compared against the wrong partner and the count explodes."""
    parent = "MKKLLAAWQTPYNARGPDPTAASLEASAGPFTVRSFTVSRP"
    construct = parent[:15] + parent[18:]          # three residues deleted internally
    assert seeds.count_substitutions(parent, construct) == 0


def test_pdb_derived_entries_are_well_formed():
    """Each PDB-derived variant names a parent that exists and a published substitution
    count to check the deposited construct against."""
    known = {w.enzyme_id for w in seeds.WILD_TYPES}
    for eid, (pdb_id, parent, expected, ref) in seeds.PDB_DERIVED.items():
        assert parent in known, f"{eid} names unknown parent {parent}"
        assert len(pdb_id) == 4, f"{eid} has a malformed PDB id {pdb_id!r}"
        assert isinstance(expected, int) and expected > 0
        assert ref
