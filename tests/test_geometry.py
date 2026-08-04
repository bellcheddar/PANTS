"""Geometry tests, including the risk-7 validation against experimental structures.

The structure-backed tests need mmCIF files fetched from RCSB into data/interim, so they
skip when those are absent rather than failing. They are the ones that matter: geometry
code that produces plausible numbers on a metagenomic prediction tells you nothing, and
only a structure with a published answer can show it is right.
"""

from __future__ import annotations

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


def test_cleft_width_orders_polyesterases_above_cutinases():
    """The discriminating feature spec section 6 predicts, tested on known structures.

    IsPETase's cleft is wider than a classic fungal cutinase's, and that ordering is the
    whole reason the structure stage can discriminate where sequence annotation cannot.
    """
    ispetase = geometry.measure(_cif("6EQE")).cleft_width_A
    lcc = geometry.measure(_cif("4EB0")).cleft_width_A
    cutinase = geometry.measure(_cif("1CEX")).cleft_width_A
    assert ispetase is not None and cutinase is not None
    assert ispetase > lcc > cutinase


def test_ispetase_aromatic_clamp_contains_the_mobile_tryptophan():
    """W185 is the residue whose mobility distinguishes IsPETase from cutinases whose
    equivalent position is fixed (spec section 2, point 1)."""
    site = geometry.measure(_cif("6EQE"))
    assert "TRP185" in site.aromatic_clamp
    assert "TYR87" in site.aromatic_clamp
