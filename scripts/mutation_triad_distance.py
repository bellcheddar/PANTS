"""How far each substituted residue sits from the catalytic triad.

The question the structural overlay grid exists to ask: do stabilising mutations keep clear
of the catalytic machinery, or do some crowd it? That is a distance, and distances belong
in a batch job rather than in a request handler -- the web process is deliberately free of
the scientific stack, and this is read from coordinates already on disk.

Distance is measured SIDE CHAIN to SIDE CHAIN, and that detail is the whole measurement.

A first pass used any heavy atom to any heavy atom and produced 1.31 to 1.35 A for four
variants -- impossible for a contact, because that is a peptide bond length. Those
mutations simply sit NEXT to a triad residue in sequence, so the minimum was picking up
backbone connectivity: W159 is bonded to catalytic S160, and a covalent bond is not
evidence that a mutation crowds the active site.

Backbone atoms (N, CA, C, O) are therefore excluded on both sides, leaving the functional
groups whose contact would actually mean something. Glycine and alanine, which have little
or no side chain, fall back to CB or CA and are flagged. Sequence separation to the nearest
triad residue is stored alongside, so an adjacent-in-sequence mutation is visible as such
rather than hiding behind a small number.

Numbering: mutation positions are stated against the lineage wild type in sequence
numbering, while the coordinates carry the depositor's. seq_offset reconciles them, and a
structure without one is skipped rather than measured in the wrong frame.
"""

from __future__ import annotations

import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config
from pipeline.db import connect, retry_write
from pipeline.structure import reference

REF_DIR = reference.REF_DIR


def distances(path: pathlib.Path, triad: list, mut_resnums: list) -> dict:
    import gemmi
    st = gemmi.read_structure(str(path))
    st.setup_entities()
    st.remove_ligands_and_waters()
    chain = st[0][0]
    by_num = {r.seqid.num: r for r in chain}

    BACKBONE = {"N", "CA", "C", "O", "OXT"}

    def side_chain(res):
        """Heavy side-chain atoms, falling back to CB then CA for the small residues."""
        atoms = [a for a in res if a.element.name != "H" and a.name not in BACKBONE]
        if atoms:
            return atoms, False
        fallback = [a for a in res if a.name in ("CB", "CA")]
        return fallback, True

    triad_atoms = []
    for t in triad:
        r = by_num.get(t)
        if r:
            triad_atoms += [(a.pos, t) for a in side_chain(r)[0]]
    if not triad_atoms:
        return {}

    out = {}
    for m in mut_resnums:
        r = by_num.get(m)
        if r is None:
            continue
        atoms, approx = side_chain(r)
        best = None
        for a in atoms:
            for pos, tnum in triad_atoms:
                d = a.pos.dist(pos)
                if best is None or d < best[0]:
                    best = (d, tnum)
        if best is None:
            continue
        sep = min(abs(m - t) for t in triad)
        out[str(m)] = {"d": round(best[0], 2), "res": r.name, "near": best[1],
                       "sep": sep, "approx": approx}
    return out


if __name__ == "__main__":
    from pipeline.recall import seeds
    muts = {v.enzyme_id: v.mutations for v in seeds.VARIANTS if v.mutations}

    with connect() as c:
        rows = c.execute(
            "SELECT rs.enzyme_id, rs.coord_path, rs.seq_offset, rg.triad_ser_resnum, "
            "       rg.triad_his_resnum, rg.triad_asp_resnum "
            "FROM reference_structures rs "
            "LEFT JOIN reference_geometry rg ON rg.enzyme_id=rs.enzyme_id "
            "WHERE rs.coord_path IS NOT NULL").fetchall()

    done = skipped = 0
    for eid, coord, off, s, h, a in rows:
        if eid not in muts or s is None:
            continue
        if off is None:
            print(f"  SKIP {eid}: no sequence offset, cannot place mutations in the structure")
            skipped += 1
            continue
        positions = [int("".join(ch for ch in m if ch.isdigit())) + off
                     for m in muts[eid] if any(ch.isdigit() for ch in m)]
        d = distances(REF_DIR / coord, [s, h, a], positions)
        if not d:
            skipped += 1
            continue

        def _do(eid=eid, d=d):
            with connect() as c:
                c.execute("UPDATE reference_structures SET mutation_geometry_json=? "
                          "WHERE enzyme_id=?", (json.dumps(d), eid))
        retry_write(_do)
        near = sorted(d.values(), key=lambda x: x["d"])[:2]
        bits = ", ".join(f"{n['res']} {n['d']} A (Δseq {n['sep']})" for n in near)
        print(f"  {eid:<32}{len(d):>3} mutations  closest {bits}")
        done += 1
    print(f"\n  measured {done}, skipped {skipped}")
