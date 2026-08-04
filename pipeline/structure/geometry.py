"""Active-site geometry: the measurement that does not inherit an annotation.

This is why the structure stage matters more than its position in the plan suggests. Every
sequence-derived label in this project traces back, directly or by similarity, to somebody
else's annotation: the Phase 5 experiment showed a head trained on those labels scores
AUC 1.000 by reproducing a similarity rule with a similarity representation. Cleft width
and the aromatic clamp are different in kind. They are measured off coordinates, so they
are the first signal here that is genuinely independent of the annotation.

Spec section 6.3 asks for triad distances and angles, oxyanion hole geometry, cleft width
and depth, aromatic residues lining the cleft, and solvent accessibility.

Validated against experimental structures where the answer is known BEFORE being run on
any prediction (PLAN_v1.md risk 7): IsPETase's catalytic triad has a published hydrogen
bonding network, so if this code cannot recover it from 6EQE it is wrong, and no amount of
plausible-looking output on a metagenomic candidate would reveal that.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

# Hydrogen bond heavy-atom distance. A catalytic triad whose Ser OG to His NE2 distance
# exceeds this is not a functional charge relay, whatever the sequence says.
HBOND_MAX_A = 3.5

# Residues forming the aromatic clamp that distinguishes polyesterases from
# soluble-ester-only esterases (spec section 6).
AROMATIC = {"TRP", "TYR", "PHE"}

# Radius around the catalytic serine within which cleft-lining residues are counted.
#
# 16 A, not 12, and the width below is a PERCENTILE rather than a maximum. Both changes
# come from a specific failure: measuring IsPETase's crystal (6EQE) gave a 20.90 A cleft
# and its own ESMFold prediction gave 11.09 A, which would have placed the reference
# PETase in cutinase territory. The fold was fine. PHE201 sits 11.7 A from the nucleophile
# in the crystal and 12.1 A in the prediction, so a 0.4 A coordinate shift moved one
# residue across a hard cutoff, and because the width was a max over pairwise distances
# that single membership change halved the answer.
#
# A hard cutoff feeding a max is brittle by construction: it makes the metric depend on a
# knife-edge decision about one residue. Both numbers looked plausible, so nothing would
# have flagged it downstream.
CLEFT_RADIUS_A = 16.0

# Percentile of pairwise aromatic separations used as the width. The max is what made the
# metric hostage to a single residue; a high percentile keeps the same meaning (how far
# apart the clamp residues sit) without one outlier deciding it.
CLEFT_WIDTH_PERCENTILE = 90


@dataclass
class ActiveSite:
    ser_resnum: Optional[int] = None
    asp_resnum: Optional[int] = None
    his_resnum: Optional[int] = None
    ser_og_his_ne2_A: Optional[float] = None
    his_nd1_asp_od_A: Optional[float] = None
    ser_his_asp_angle_deg: Optional[float] = None
    oxyanion_n1_A: Optional[float] = None
    oxyanion_n2_A: Optional[float] = None
    cleft_width_A: Optional[float] = None
    cleft_depth_A: Optional[float] = None
    aromatic_clamp: List[str] = field(default_factory=list)
    n_cleft_residues: Optional[int] = None
    n_aromatic_lining: Optional[int] = None
    # Kept alongside the percentile so the brittle quantity stays visible: a large gap
    # between the two means one residue is dominating and the measure is unstable there.
    cleft_width_max_A: Optional[float] = None

    @property
    def triad_is_connected(self) -> bool:
        """Both charge-relay hydrogen bonds present. A triad in sequence is not a triad
        in space unless the distances say so."""
        return (self.ser_og_his_ne2_A is not None and self.ser_og_his_ne2_A <= HBOND_MAX_A
                and self.his_nd1_asp_od_A is not None and self.his_nd1_asp_od_A <= HBOND_MAX_A)


def load_structure(path: str | Path):
    """Read a PDB or mmCIF into a gemmi structure, first model only."""
    import gemmi
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_ligands_and_waters()
    return st


def _atom(residue, *names):
    for n in names:
        a = residue.find_atom(n, "*")
        if a is not None:
            return a
    return None


def _pos(atom) -> Optional[np.ndarray]:
    return None if atom is None else np.array([atom.pos.x, atom.pos.y, atom.pos.z])


def _dist(a, b) -> Optional[float]:
    if a is None or b is None:
        return None
    return float(np.linalg.norm(a - b))


def _angle(a, b, c) -> Optional[float]:
    if a is None or b is None or c is None:
        return None
    v1, v2 = a - b, c - b
    cosang = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2) + 1e-12))
    return float(math.degrees(math.acos(max(-1.0, min(1.0, cosang)))))


def find_triad(structure, chain_id: Optional[str] = None) -> Optional[Tuple]:
    """Locate the Ser-His-Asp triad by GEOMETRY, not by sequence position.

    Searching for a spatially connected Ser/His/Asp is what makes this independent of the
    sequence-level annotation: a candidate whose three residues are present but not in
    contact does not have a catalytic triad, and only coordinates can say so.
    """
    chain = None
    for ch in structure[0]:
        if chain_id is None or ch.name == chain_id:
            chain = ch
            break
    if chain is None:
        return None

    sers, hiss, asps = [], [], []
    for res in chain:
        if res.name == "SER":
            sers.append(res)
        elif res.name == "HIS":
            hiss.append(res)
        elif res.name in ("ASP", "GLU"):
            asps.append(res)

    best = None
    for ser in sers:
        og = _pos(_atom(ser, "OG"))
        if og is None:
            continue
        for his in hiss:
            ne2 = _pos(_atom(his, "NE2"))
            nd1 = _pos(_atom(his, "ND1"))
            if ne2 is None or nd1 is None:
                continue
            d1 = _dist(og, ne2)
            if d1 is None or d1 > HBOND_MAX_A:
                continue
            for asp in asps:
                od = _pos(_atom(asp, "OD1", "OD2", "OE1", "OE2"))
                d2 = _dist(nd1, od)
                if d2 is None or d2 > HBOND_MAX_A:
                    continue
                score = d1 + d2
                if best is None or score < best[0]:
                    best = (score, ser, his, asp, d1, d2)
    return best


def measure(path: str | Path, chain_id: Optional[str] = None) -> ActiveSite:
    """Full active-site measurement for one structure."""
    st = load_structure(path)
    site = ActiveSite()
    found = find_triad(st, chain_id)
    if found is None:
        return site

    _score, ser, his, asp, d1, d2 = found
    site.ser_resnum = ser.seqid.num
    site.his_resnum = his.seqid.num
    site.asp_resnum = asp.seqid.num
    site.ser_og_his_ne2_A = round(d1, 2)
    site.his_nd1_asp_od_A = round(d2, 2)
    site.ser_his_asp_angle_deg = None

    og = _pos(_atom(ser, "OG"))
    ne2 = _pos(_atom(his, "NE2"))
    od = _pos(_atom(asp, "OD1", "OD2", "OE1", "OE2"))
    a = _angle(og, ne2, od)
    if a is not None:
        site.ser_his_asp_angle_deg = round(a, 1)

    chain = next((ch for ch in st[0] if chain_id is None or ch.name == chain_id), None)
    if chain is None:
        return site
    residues = list(chain)

    # --- oxyanion hole: backbone amide nitrogens nearest the nucleophile ---
    # Structural, not sequence-derived. The two closest backbone N atoms (excluding the
    # serine's own) are the donors that stabilise the transition state.
    n_dists: List[Tuple[float, int]] = []
    for res in residues:
        if res.seqid.num == ser.seqid.num:
            continue
        n = _pos(_atom(res, "N"))
        d = _dist(og, n)
        if d is not None and d < 8.0:
            n_dists.append((d, res.seqid.num))
    n_dists.sort()
    if n_dists:
        site.oxyanion_n1_A = round(n_dists[0][0], 2)
    if len(n_dists) > 1:
        site.oxyanion_n2_A = round(n_dists[1][0], 2)

    # --- cleft: residues lining the pocket around the nucleophile ---
    lining: List[Tuple[float, object]] = []
    for res in residues:
        ca = _pos(_atom(res, "CA"))
        d = _dist(og, ca)
        if d is not None and d <= CLEFT_RADIUS_A:
            lining.append((d, res))
    site.n_cleft_residues = len(lining)
    site.aromatic_clamp = sorted(
        f"{r.name}{r.seqid.num}" for d, r in lining if r.name in AROMATIC)

    # Cleft width: the widest separation among aromatic side-chain centroids lining the
    # pocket. IsPETase's cleft is wider than a cutinase's, which is the discriminating
    # feature spec section 6 asks to surface.
    centroids = []
    for d, r in lining:
        if r.name not in AROMATIC:
            continue
        pts = [_pos(a) for a in r if a.name not in ("N", "CA", "C", "O")]
        pts = [p for p in pts if p is not None]
        if pts:
            centroids.append(np.mean(pts, axis=0))
    site.n_aromatic_lining = len(centroids)
    if len(centroids) >= 2:
        M = np.array(centroids)
        dists = np.linalg.norm(M[:, None, :] - M[None, :, :], axis=-1)
        upper = dists[np.triu_indices(len(M), k=1)]
        site.cleft_width_A = round(float(np.percentile(upper, CLEFT_WIDTH_PERCENTILE)), 2)
        site.cleft_width_max_A = round(float(upper.max()), 2)

    # Cleft depth: how far the nucleophile sits below the local surface, approximated as
    # the distance from OG to the centroid of the lining CA atoms.
    cas = [_pos(_atom(r, "CA")) for d, r in lining]
    cas = [c for c in cas if c is not None]
    if cas:
        site.cleft_depth_A = round(float(np.linalg.norm(og - np.mean(cas, axis=0))), 2)

    return site
