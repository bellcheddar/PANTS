"""Tests for profile-column mapping in triad detection.

The bug these exist to prevent: hmmalign writes residues that are insertions relative to
the model in LOWER case, and match/delete states in upper case. Both are real residues of
the sequence and both consume a position in its own numbering. An earlier version counted
only upper-case characters when locating the reference's Nth residue, which walked off by
the number of insertions before the target and made IsPETase report its own verified
catalytic Ser160 as A171.

That failure is silent in production: the triad simply reads as incomplete and the
candidate is discarded, so a broken mapping looks like a strict filter rather than a bug.
The only reason it was caught is that the reference's answer was known in advance.
"""

from __future__ import annotations

import pytest

from pipeline.recall import triad
from pipeline.recall.triad import _read_stockholm, _residue_at_model_position


def test_lowercase_insertions_count_towards_reference_numbering():
    """A lower-case insertion before the target must shift the column, not be skipped."""
    #             col: 0123456
    ref_aligned = "ACdeFGH"     # ungapped reference is A C d e F G H -> 7 residues
    # The 5th residue of the reference is 'F' (A=1, C=2, d=3, e=4, F=5), at column 4.
    res, pos = _residue_at_model_position(ref_aligned, ref_aligned, 5)
    assert res == "F"
    assert pos == 5


def test_counting_only_uppercase_would_pick_the_wrong_column():
    """Guards the exact regression: upper-case-only counting picks residue 5 as 'H'."""
    ref_aligned = "ACdeFGH"
    uppercase_only = [i for i, c in enumerate(ref_aligned) if c.isupper()]
    # Upper-case positions are A,C,F,G,H -> the 5th would be 'H' at column 6, not 'F'.
    assert ref_aligned[uppercase_only[4]] == "H"
    assert _residue_at_model_position(ref_aligned, ref_aligned, 5)[0] == "F"


def test_gap_in_candidate_reports_none():
    ref_aligned = "ACDEFGH"
    candidate = "AC-EFGH"
    res, pos = _residue_at_model_position(candidate, ref_aligned, 3)
    assert res is None and pos is None


def test_candidate_position_is_its_own_ungapped_numbering():
    """The returned position must index the CANDIDATE, not the alignment column."""
    ref_aligned = "ACDEFGH"
    candidate = "--DEFGH"      # candidate's D is its residue 1
    res, pos = _residue_at_model_position(candidate, ref_aligned, 3)
    assert res == "D"
    assert pos == 1


def test_reference_position_beyond_sequence_returns_none():
    assert _residue_at_model_position("ACDEF", "ACDEF", 99) == (None, None)


def test_triad_call_completeness():
    call = triad.TriadCall(sequence_id="x", ser="S", asp="D", his="H",
                           ser_pos=1, asp_pos=2, his_pos=3,
                           oxyanion_residues={}, aligned_fraction=1.0)
    assert call.complete
    assert call.reason == "complete"

    bad = triad.TriadCall(sequence_id="y", ser="S", asp="K", his="H",
                          ser_pos=1, asp_pos=2, his_pos=3,
                          oxyanion_residues={}, aligned_fraction=1.0)
    assert not bad.complete
    assert "Asp" in bad.reason


def test_stockholm_concatenates_interleaved_blocks(tmp_path):
    """Long alignments come back interleaved; a later block must extend, not replace."""
    sto = tmp_path / "a.sto"
    sto.write_text(
        "# STOCKHOLM 1.0\n"
        "seq1 ACDE\n"
        "seq2 WXYZ\n"
        "\n"
        "seq1 FGHI\n"
        "seq2 KLMN\n"
        "//\n"
    )
    parsed = _read_stockholm(sto)
    assert parsed["seq1"] == "ACDEFGHI"
    assert parsed["seq2"] == "WXYZKLMN"


def test_ispetase_reference_triad_is_the_verified_one():
    """These three numbers were checked against the real UniProt sequence in Phase 1.
    Changing them silently re-anchors every triad call in the project."""
    assert triad.ISPETASE_TRIAD == {"ser": 160, "asp": 206, "his": 237}


# --- anchors: the per-cluster replacement for a single hardcoded reference -------------

def test_anchor_consistency_requires_sdh():
    from pipeline.recall.anchors import Anchor
    seq = "AAASAAADAAAH"          # S at 4, D at 8, H at 12
    good = Anchor("x", "P1", seq, ser=4, asp=8, his=12, source="uniprot")
    assert good.is_consistent
    assert good.residues() == "SDH"

    bad = Anchor("y", "P2", seq, ser=1, asp=8, his=12, source="uniprot")
    assert not bad.is_consistent


def test_anchor_assigns_triad_by_residue_not_by_description():
    """UniProt's charge-relay entries do not say which member is Asp and which is His, so
    assignment must come from the residue itself."""
    import pipeline.recall.anchors as A

    seq = "MMMSMMMMDMMMH"        # S4, D9, H13
    def fake(_acc):
        # deliberately unordered, with vague descriptions
        return seq, [(13, "charge relay system"), (4, "nucleophile"), (9, "charge relay system")]

    orig, A.fetch_active_sites = A.fetch_active_sites, fake
    try:
        a = A.anchor_from_uniprot("z", "P0")
    finally:
        A.fetch_active_sites = orig

    assert a is not None
    assert (a.ser, a.asp, a.his) == (4, 9, 13)
    assert a.is_consistent


def test_anchor_rejected_when_fewer_than_three_sites():
    import pipeline.recall.anchors as A
    orig, A.fetch_active_sites = A.fetch_active_sites, lambda _a: ("MSMDMH", [(2, "nucleophile")])
    try:
        assert A.anchor_from_uniprot("z", "P0") is None
    finally:
        A.fetch_active_sites = orig
