#!/usr/bin/env python3
"""Pool the language model over the active site only, instead of the whole protein.

Every evaluation in this project has said the same thing in a different way: the labels
track family membership, and a model given the whole sequence learns the family. Mean
pooling makes that easy -- most of a 300-residue embedding describes a fold that both
classes share, so the few positions where a PET-degrader differs from a polyesterase that
cannot degrade PET are averaged into noise.

This asks the narrower question the project actually cares about. Take the residues that
line the catalytic site in the folded structure, and pool the per-residue embeddings over
those positions alone. The global fold signal is discarded by construction, and what
remains is the local sequence environment of the chemistry.

That is only possible now because the 205-structure fold finished: 434 of the 447 enzymes
in the contrast, 97%, have a structure with a measured triad and a reconciled sequence
offset. The offset matters more than it sounds. A deposit is numbered by its depositor, so
residue 160 in the coordinates is not residue 160 in the stored sequence -- three of the
first twenty references were out by +42, +26 and -40 -- and pooling over unreconciled
positions would silently average the wrong residues while producing a perfectly plausible
vector.

Two radii are tried, 8 and 12 A from the catalytic serine's OG, because "the active site"
has no single definition and a result that only appears at one radius is a result about
the radius.
"""

from __future__ import annotations

import json
import pathlib
import sys
from typing import Dict, List, Optional, Set, Tuple

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config, seqtools
from pipeline.db import connect

import importlib.util
_spec = importlib.util.spec_from_file_location(
    "shv", pathlib.Path(__file__).resolve().parent / "sequence_head_variants.py")
shv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(shv)

RADII = (8.0, 12.0)
MODEL = "facebook/esm2_t33_650M_UR50D"
REF_DIR = config.STATIC_DIR / "reference_structures"
OUT = config.INTERIM_DIR / "active_site_embeddings.json"


def _atoms(path: pathlib.Path):
    """(resnum, atom name, xyz) for the first altloc of every ATOM record."""
    out = []
    for line in path.read_text().splitlines():
        if line.startswith("ATOM") and line[16] in " A":
            try:
                out.append((int(line[22:26]), line[12:16].strip(),
                            (float(line[30:38]), float(line[38:46]), float(line[46:54]))))
            except ValueError:
                continue
    return out


def site_positions(enzyme_id: str, coord_path: str, ser_resnum: int, offset: int,
                   radius: float, seq_len: int) -> Set[int]:
    """1-based SEQUENCE positions whose residue lines the active site.

    Centred on the catalytic serine's OG where present, its CA otherwise, and returned in
    sequence numbering: `structure_resnum = sequence_position + offset`, so the inverse is
    applied here. Positions outside the stored sequence are dropped rather than clamped --
    a deposit can model residues the stored sequence does not contain.
    """
    path = REF_DIR / coord_path
    if not path.exists():
        return set()
    atoms = _atoms(path)
    centre = None
    for resnum, name, xyz in atoms:
        if resnum == ser_resnum and name == "OG":
            centre = xyz
            break
    if centre is None:
        for resnum, name, xyz in atoms:
            if resnum == ser_resnum and name == "CA":
                centre = xyz
                break
    if centre is None:
        return set()
    cx, cy, cz = centre
    r2 = radius * radius
    near = {resnum for resnum, name, (x, y, z) in atoms
            if (x - cx) ** 2 + (y - cy) ** 2 + (z - cz) ** 2 <= r2}
    return {p for p in (n - offset for n in near) if 1 <= p <= seq_len}


def load_sites(radius: float) -> Dict[str, Set[int]]:
    with connect() as c:
        rows = c.execute(
            "SELECT rs.enzyme_id, rs.coord_path, rs.seq_offset, rg.triad_ser_resnum, "
            "       ce.seq_length "
            "FROM reference_structures rs JOIN reference_geometry rg USING(enzyme_id) "
            "JOIN characterised_enzymes ce USING(enzyme_id) "
            "WHERE rg.triad_ser_resnum IS NOT NULL AND rs.coord_path IS NOT NULL "
            "  AND rs.seq_offset IS NOT NULL").fetchall()
    out = {}
    for eid, coord, offset, ser, seq_len in rows:
        pos = site_positions(eid, coord, ser, offset or 0, radius, seq_len or 0)
        if pos:
            out[eid] = pos
    return out


def embed_pooled(recs, sites: Dict[str, Set[int]], model_name: str) -> Tuple[np.ndarray, List[str]]:
    """Mean over the site positions of each sequence's per-residue embeddings."""
    import torch
    from pipeline.embed import esm

    tok, model = esm.load_model(model_name)
    vecs, kept = [], []
    for n, (eid, seq) in enumerate(recs):
        pos = sites.get(eid)
        if not pos:
            continue
        enc = tok(seq[:esm.MAX_LENGTH], return_tensors="pt", truncation=True,
                  max_length=esm.MAX_LENGTH)
        out = model(**enc).last_hidden_state[0]
        # Token 0 is CLS, so sequence position p is token p. Anything past truncation is
        # simply absent rather than mapped to the wrong token.
        idx = [p for p in sorted(pos) if p < out.shape[0]]
        if not idx:
            continue
        vecs.append(out[idx].mean(0).numpy())
        kept.append(eid)
        if n % 100 == 0:
            print(f"    pooled {n}/{len(recs)}", flush=True)
    return np.vstack(vecs), kept


def main() -> int:
    recs, y = shv.contrast()
    by_id = dict(recs)
    label = {e: int(v) for (e, _), v in zip(recs, y)}
    results = {}

    for radius in RADII:
        sites = load_sites(radius)
        n_res = [len(v) for v in sites.values()]
        print(f"\nradius {radius:g} A: {len(sites)} enzymes, "
              f"median {int(np.median(n_res))} residues in the site")

        sub = [(e, s) for e, s in recs if e in sites]
        X, kept = embed_pooled(sub, sites, MODEL)
        yk = np.array([label[e] for e in kept])
        keep, groups = shv.mixed_cluster_mask([(e, by_id[e]) for e in kept], yk)
        Xm, ym, gm = X[keep], yk[keep], groups[keep]
        if len(set(ym)) < 2:
            print("    not evaluable: one class after restriction")
            continue
        comp = shv.composition([(e, by_id[e]) for e, k in zip(kept, keep) if k])
        c_mean, c_sd = shv.grouped_auc(comp, ym, gm)
        mean, sd = shv.grouped_auc(Xm, ym, gm)
        print(f"    mixed clusters: {len(ym)} enzymes, {int(ym.sum())} active, "
              f"{int((1-ym).sum())} inactive, {len(set(gm))} clusters")
        print(f"    active-site pooled : {mean:.3f} +/- {sd:.3f}")
        print(f"    composition        : {c_mean:.3f} +/- {c_sd:.3f}")
        results[f"radius_{radius:g}A"] = {
            "auc": mean, "sd": sd, "composition": c_mean, "n": len(ym),
            "n_clusters": len(set(gm)), "median_site_residues": int(np.median(n_res))}

    OUT.write_text(json.dumps(results, indent=2))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
