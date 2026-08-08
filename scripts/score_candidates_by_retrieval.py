#!/usr/bin/env python3
"""Rank candidates by similarity to the nearest MEASURED-ACTIVE enzyme, and say when that
ranking is worth anything.

This replaces a learned classifier as the ranking method, on evidence. Cross-lineage, a
trained head is at chance and its coefficients point in opposite directions between
lineages; within a lineage it cannot be shown to beat nearest-neighbour retrieval
(paired lead +0.159 [-0.263, +0.390]). Retrieval needs no training, cannot learn a lineage
direction that reverses on an unseen family, and nothing built here has beaten it on labels
somebody measured.

Two things separate this from the `nearest_characterised_id` the recall stage already
records. That field points at the nearest enzyme in the WHOLE catalogue, most of which
carry EC numbers assigned automatically by similarity -- so it can rank a candidate highly
for resembling something nobody has ever assayed. This one restricts the reference set to
enzymes whose PET activity was MEASURED and published.

And it reports its own competence. The identity-decay curve says the signal is weak at 70%
identity and unmeasurable below 50%, so a score is only meaningful for a candidate close
enough to a measured enzyme for that to apply. Every candidate is therefore banded:

    in-range   >= 70% identity to a measured-active enzyme
    marginal   50-70%
    out-of-range < 50%   -- ranked, but the ranking has no demonstrated validity here

Publishing the band alongside the score is the point. A ranked list with no statement of
where it stops working is the failure mode this whole project has been arguing against.
"""

from __future__ import annotations

import json
import pathlib
import sys
from typing import Dict, Tuple

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config, seqtools
from pipeline.db import connect, now

MODEL_VERSION = "retrieval/nearest-measured-active/v1"
IN_RANGE, MARGINAL = 0.70, 0.50
OUT = config.INTERIM_DIR / "candidate_retrieval_scores.json"


def band(identity: float) -> str:
    if identity >= IN_RANGE:
        return "in-range"
    if identity >= MARGINAL:
        return "marginal"
    return "out-of-range"


def main() -> int:
    with connect() as c:
        cands = c.execute(
            "SELECT candidate_id, sequence FROM candidates WHERE sequence IS NOT NULL").fetchall()
        marks = ",".join("?" * len(config.MEASURED_TIERS))
        refs = c.execute(
            f"SELECT enzyme_id, sequence FROM characterised_enzymes "
            f"WHERE is_positive=1 AND source_ref IN ({marks}) AND excluded_from_training=0 "
            f"  AND sequence IS NOT NULL AND family!='mhetase_like'",
            config.MEASURED_TIERS).fetchall()
    print(f"{len(cands)} candidates against {len(refs)} MEASURED-ACTIVE reference enzymes")

    q = seqtools.write_fasta([tuple(r) for r in cands], config.INTERIM_DIR / "cand_q.fasta")
    t = seqtools.write_fasta([tuple(r) for r in refs], config.INTERIM_DIR / "ref_measured.fasta")
    best: Dict[str, Tuple[str, float]] = seqtools.best_identity_to(q, t)
    print(f"{len(best)} candidates matched a reference at all\n")

    rows, counts = [], {"in-range": 0, "marginal": 0, "out-of-range": 0, "no match": 0}
    for cid, _ in cands:
        hit = best.get(cid)
        if not hit:
            counts["no match"] += 1
            rows.append({"candidate_id": cid, "nearest": None, "identity": None,
                         "band": "no match"})
            continue
        ref, ident = hit
        b = band(ident)
        counts[b] += 1
        rows.append({"candidate_id": cid, "nearest": ref, "identity": ident, "band": b})

    with connect() as c:
        c.executemany(
            "INSERT INTO scores(candidate_id, pet_activity_prob, model_version, scored_at) "
            "VALUES (?,?,?,?) ON CONFLICT(candidate_id) DO UPDATE SET "
            " pet_activity_prob=excluded.pet_activity_prob, "
            " model_version=excluded.model_version, scored_at=excluded.scored_at",
            [(r["candidate_id"], r["identity"], MODEL_VERSION, now())
             for r in rows if r["identity"] is not None])
        # The nearest-measured-active neighbour, kept beside the catalogue's existing
        # nearest-ANY-enzyme so the two are not confused for each other.
        c.executemany(
            "UPDATE candidates SET nearest_measured_id=?, nearest_measured_identity=?, "
            " retrieval_band=? WHERE candidate_id=?",
            [(r["nearest"], r["identity"], r["band"], r["candidate_id"]) for r in rows])
        c.commit()

    print(f"{'band':<14} {'candidates':>10}   meaning")
    print(f"{'in-range':<14} {counts['in-range']:>10}   >= 70% identity to a measured-active enzyme")
    print(f"{'marginal':<14} {counts['marginal']:>10}   50-70%")
    print(f"{'out-of-range':<14} {counts['out-of-range']:>10}   < 50%, ranking has no demonstrated validity")
    print(f"{'no match':<14} {counts['no match']:>10}   no alignment to any measured-active enzyme")

    idents = sorted(r["identity"] for r in rows if r["identity"] is not None)
    if idents:
        import statistics
        print(f"\nidentity to nearest measured-active enzyme: "
              f"median {statistics.median(idents):.1%}, max {max(idents):.1%}")

    OUT.write_text(json.dumps({"model_version": MODEL_VERSION, "counts": counts,
                               "n_references": len(refs), "rows": rows}, indent=2))
    print(f"\nwrote {OUT}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
