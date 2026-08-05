"""Geometry tests, including the risk-7 validation against experimental structures.

The structure-backed tests need mmCIF files fetched from RCSB into data/interim, so they
skip when those are absent rather than failing. They are the ones that matter: geometry
code that produces plausible numbers on a metagenomic prediction tells you nothing, and
only a structure with a published answer can show it is right.
"""

from __future__ import annotations

import pathlib

import numpy as np
import pytest

from pipeline import config
from pipeline.structure import geometry


def _cif(pdb_id: str):
    p = config.INTERIM_DIR / "pdb" / f"{pdb_id}.cif"
    if not p.exists():
        pytest.skip(f"{pdb_id}.cif not fetched")
    return p


# --- pure geometry helpers -------------------------------------------------------------

def test_distance_and_angle():
    a = np.array([0.0, 0.0, 0.0])
    b = np.array([1.0, 0.0, 0.0])
    c = np.array([1.0, 1.0, 0.0])
    assert geometry._dist(a, b) == pytest.approx(1.0)
    assert geometry._angle(a, b, c) == pytest.approx(90.0)


def test_distance_with_missing_atom_is_none():
    assert geometry._dist(None, np.zeros(3)) is None
    assert geometry._angle(np.zeros(3), None, np.zeros(3)) is None


def test_triad_connectivity_requires_both_hydrogen_bonds():
    """Three residues present in sequence is not a triad. The distances decide."""
    connected = geometry.ActiveSite(ser_og_his_ne2_A=2.9, his_nd1_asp_od_A=3.0)
    assert connected.triad_is_connected

    too_far = geometry.ActiveSite(ser_og_his_ne2_A=2.9, his_nd1_asp_od_A=6.4)
    assert not too_far.triad_is_connected

    missing = geometry.ActiveSite(ser_og_his_ne2_A=2.9, his_nd1_asp_od_A=None)
    assert not missing.triad_is_connected


# --- risk 7: known answers from experimental structures --------------------------------

@pytest.mark.parametrize("pdb_id,ser,asp,his", [
    ("6EQE", 160, 206, 237),   # IsPETase, UniProt precursor numbering
    ("4EB0", 165, 210, 242),   # LCC
])
def test_recovers_published_triad(pdb_id, ser, asp, his):
    site = geometry.measure(_cif(pdb_id))
    assert (site.ser_resnum, site.asp_resnum, site.his_resnum) == (ser, asp, his)
    assert site.triad_is_connected


def test_tfcut2_offset_is_the_signal_peptide_not_an_error():
    """4CG1 numbers the MATURE protein; UniProt numbers the precursor.

    The triad reads S130/D176/H208 against UniProt's S170/D216/H248: a constant offset of
    40, which is exactly TfCut2's signal peptide length (1-40). A constant offset across
    all three residues is the signature of a numbering convention, not a detection error,
    and this test exists so nobody 'fixes' it later.
    """
    site = geometry.measure(_cif("4CG1"))
    assert (site.ser_resnum, site.asp_resnum, site.his_resnum) == (130, 176, 208)
    offsets = {170 - site.ser_resnum, 216 - site.asp_resnum, 248 - site.his_resnum}
    assert offsets == {40}


def test_cleft_width_separates_plc_family_from_fungal_cutinase():
    """What the robust metric actually supports, which is less than first claimed.

    An earlier version used a MAX over pairwise aromatic separations inside a hard 12 A
    radius, and reported IsPETase 20.9 > LCC 17.8 > TfCut2 12.8 > cutinase 11.2: a neat
    two-fold range in the expected order. That ordering was an artefact. The max made the
    measure hostage to a single residue crossing the radius, and IsPETase's own ESMFold
    prediction scored 11.09 against its crystal's 20.90 for exactly that reason.

    With the percentile measure the honest claim is narrower: the classic FUNGAL cutinase
    (1CEX) is clearly narrower than the Polyesterase-lipase-cutinase family members, but
    those members do not order among themselves by PET activity. Separating families is
    not the same as separating PET-active from PET-inactive within a family, and only the
    former is demonstrated here.
    """
    ispetase = geometry.measure(_cif("6EQE")).cleft_width_A
    lcc = geometry.measure(_cif("4EB0")).cleft_width_A
    tfcut2 = geometry.measure(_cif("4CG1")).cleft_width_A
    cutinase = geometry.measure(_cif("1CEX")).cleft_width_A
    assert min(ispetase, lcc, tfcut2) > cutinase + 3.0


def test_cleft_width_is_stable_between_crystal_and_prediction():
    """The regression that motivated the percentile: a 0.4 A shift in one residue must not
    move the answer by 9 A."""
    pred = config.INTERIM_DIR / "structures" / "IsPETase_esmfold.pdb"
    if not pred.exists():
        pytest.skip("ESMFold prediction not present")
    xtal = geometry.measure(_cif("6EQE")).cleft_width_A
    predicted = geometry.measure(pred).cleft_width_A
    assert abs(xtal - predicted) < 4.0, (
        f"crystal {xtal} vs prediction {predicted}: the measure is unstable again")


def test_ispetase_aromatic_clamp_contains_the_mobile_tryptophan():
    """W185 is the residue whose mobility distinguishes IsPETase from cutinases whose
    equivalent position is fixed (spec section 2, point 1)."""
    site = geometry.measure(_cif("6EQE"))
    assert "TRP185" in site.aromatic_clamp
    assert "TYR87" in site.aromatic_clamp


# --- recall robustness -----------------------------------------------------------------

def test_hmmscan_returns_nothing_for_an_empty_candidate_set():
    """An assembly whose prefilter finds nothing is a normal outcome, not an error.

    hmmscan exits non-zero on an empty input file ("empty or misformatted"), which crashed
    a 50-file gut scan on file 46. Several gut assemblies genuinely contain no
    polyesterase-like protein, so this path is common rather than exotic.
    """
    from pipeline.recall import library
    lib = library.Library(entries={}, db_path=None)
    assert library.hmmscan_best(lib, []) == {}
    assert library.call_triads(lib, []) == ({}, {})


# --------------------------------------------------------------------------------------
# Oxyanion hole: validated against a structure whose answer is published
# --------------------------------------------------------------------------------------

REFERENCE = pathlib.Path(__file__).resolve().parent.parent / "app/static/reference/6EQE.cif"


@pytest.mark.skipif(not REFERENCE.exists(), reason="6EQE reference not present")
def test_oxyanion_hole_recovers_the_published_ispetase_donors():
    """IsPETase's oxyanion hole is formed by the backbone amides of Tyr87 and Met161.

    This is the check the module docstring demands: if the code cannot recover a known
    answer from a crystal structure, no amount of plausible output on a metagenomic
    prediction would reveal that it is wrong. The previous rule -- "the two backbone N
    atoms closest to OG" -- failed it, returning Met161 and Trp185, because Tyr87 is only
    the third closest.
    """
    site = geometry.measure(REFERENCE)
    assert {site.oxyanion_n1_resnum, site.oxyanion_n2_resnum} == {87, 161}


@pytest.mark.skipif(not REFERENCE.exists(), reason="6EQE reference not present")
def test_oxyanion_rejects_the_nearer_non_donor():
    """Trp185 is closer to the nucleophile than Tyr87 and is NOT a donor.

    It is excluded on direction, not distance: its backbone N-H points away from the
    pocket (76 degrees, against the 60 degree cutoff). Pinning this stops anyone
    "simplifying" the angle test back out.
    """
    site = geometry.measure(REFERENCE)
    assert site.oxyanion_n2_resnum != 185
    assert site.oxyanion_n2_angle_deg is not None
    assert site.oxyanion_n2_angle_deg <= geometry.OXYANION_MAX_ANGLE_DEG


@pytest.mark.skipif(not REFERENCE.exists(), reason="6EQE reference not present")
def test_oxyanion_donor_one_is_the_nucleophile_elbow():
    """Donor 1 is taken by position, not searched for: it is always the residue after the
    nucleophile in an alpha/beta-hydrolase."""
    site = geometry.measure(REFERENCE)
    assert site.oxyanion_n1_resnum == site.ser_resnum + 1


def test_reference_plddt_is_normalised_to_the_same_scale():
    """AlphaFold writes pLDDT 0-100 into the B-factor column and ESMFold writes 0-1.

    Stored unreconciled they landed in one database column together, so an ESMFold model
    read 0.96 beside an AlphaFold one at 92.01 and the page rendered "1.0" as a confidence
    score. The scale is detected from the values, not assumed from which source produced
    them.
    """
    from pipeline.structure import reference

    def pdb(bfactors):
        return "\n".join(
            f"ATOM  {i:5d}  CA  ALA A{i:4d}    "
            f"   0.000   0.000   0.000  1.00{b:6.2f}           C"
            for i, b in enumerate(bfactors, start=1))

    assert reference._mean_plddt(pdb([0.80, 0.90, 1.00])) == 90.0   # 0-1 input, rescaled
    assert reference._mean_plddt(pdb([80.0, 90.0, 100.0])) == 90.0  # already 0-100, left alone
    assert reference._mean_plddt("no atoms here") is None


# --------------------------------------------------------------------------------------
# Catalytic knockout detection
# --------------------------------------------------------------------------------------

def _kn(pdb_seq, ref_seq, triad):
    from pipeline.structure import reference
    return reference.catalytic_knockout(pdb_seq, ref_seq, triad)


def test_knockout_detects_a_mutated_nucleophile():
    """The pattern this project adopted three times: a deposit carrying the right name and
    the right length with the catalytic serine mutated to alanine so the enzyme could be
    crystallised holding substrate."""
    ref = "MKTAYSGGHLDPQSER"
    ko = ref.replace("S", "A", 1)          # first serine knocked out
    hits = _kn(ko, ref, {"ser": 6})        # position 6 is the S in ...AYSGG...
    assert hits and hits[0]["role"] == "ser" and hits[0]["found"] == "A"


def test_knockout_passes_an_intact_construct():
    ref = "MKTAYSGGHLDPQSER"
    assert _kn(ref, ref, {"ser": 6, "his": 9, "asp": 11}) == []


def test_knockout_accepts_glutamate_for_the_acid():
    """Asp or Glu both complete a charge relay; flagging Glu would be a false positive."""
    ref = "MKTAYSGGHLDPQSER"
    swapped = ref[:10] + "E" + ref[11:]
    assert _kn(swapped, ref, {"asp": 11}) == []


def test_knockout_tolerates_a_different_construct_numbering():
    """The reference is the lineage wild type, and a deposit may be a shorter mature
    construct. Comparing by position rather than by alignment reported HotPETase as a
    triple knockout purely because its stored form starts 26 residues later."""
    ref = "MMMMMMMMMMKTAYSGGHLDPQSER"      # precursor, triad at 15/18/20
    mature = ref[10:]                       # same protein, numbering shifted by 10
    assert _kn(mature, ref, {"ser": 15, "his": 18, "asp": 20}) == []


def test_knockout_reports_an_unmodelled_catalytic_residue():
    ref = "MKTAYSGGHLDPQSER"
    assert _kn("MKTAY", ref, {"his": 9})[0]["note"] == "not modelled"
