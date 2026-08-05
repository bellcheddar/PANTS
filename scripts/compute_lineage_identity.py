"""Percent identity of each named enzyme to the wild type at the root of its lineage.

Precomputed rather than done in the web layer: the alignment needs biotite, and the web
virtual environment is kept free of the scientific stack because the droplet has 3.9 GB
shared with five other applications.

Identity is measured against the LINEAGE wild type, not against IsPETase, and the two are
very different questions. LCC-ICCG is about 50% identical to IsPETase and 98.6% identical
to LCC, and only the second number says anything about how far engineering has moved it.
Wild types are 100% against themselves by definition, which is honest for a lineage table.

Gapped alignment, and only positions where both sequences have a residue are counted, so a
mature construct that is shorter than its precursor parent is not penalised for the
residues it never had.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline.db import connect, retry_write


def lineage_root(enzyme_id: str, parent_of: dict) -> str:
    """Walk up matched_positive_id to the first enzyme with no parent."""
    seen = {enzyme_id}
    cur = enzyme_id
    while True:
        nxt = parent_of.get(cur)
        if not nxt or nxt in seen:      # no parent, or a cycle
            return cur
        seen.add(nxt)
        cur = nxt


def identity(a: str, b: str) -> float:
    import biotite.sequence as bseq
    import biotite.sequence.align as balign
    matrix = balign.SubstitutionMatrix.std_protein_matrix()
    aln = balign.align_optimal(bseq.ProteinSequence(a), bseq.ProteinSequence(b),
                               matrix, gap_penalty=(-10, -1))[0]
    x, y = balign.get_symbols(aln)
    both = [(p, q) for p, q in zip(x, y) if p and q]
    if not both:
        return 0.0
    return 100.0 * sum(1 for p, q in both if p == q) / len(both)


if __name__ == "__main__":
    with connect() as c:
        rows = {r[0]: (r[1], r[2]) for r in c.execute(
            "SELECT enzyme_id, sequence, matched_positive_id FROM characterised_enzymes "
            "WHERE enzyme_id NOT LIKE '%:%' AND sequence IS NOT NULL")}
    parent_of = {k: v[1] for k, v in rows.items() if v[1]}

    print(f"{'enzyme':<24}{'lineage WT':<16}{'%ID':>7}")
    for eid, (seq, _p) in sorted(rows.items()):
        root = lineage_root(eid, parent_of)
        root_seq = rows.get(root, (None, None))[0]
        if not root_seq:
            continue
        pid = 100.0 if root == eid else identity(seq, root_seq)
        print(f"{eid:<24}{root:<16}{pid:>7.1f}")

        def _do(eid=eid, root=root, pid=pid):
            with connect() as c:
                c.execute("UPDATE characterised_enzymes SET lineage_wt_id=?, "
                          "identity_to_lineage_wt=? WHERE enzyme_id=?",
                          (root, round(pid, 2), eid))
        retry_write(_do)
