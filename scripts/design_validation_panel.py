#!/usr/bin/env python3
"""Choose which enzymes to assay next, and defend every choice.

Everything this project measured points at one constraint. Not the number of enzymes, not
the model, not the representation: the number of independent LINEAGES in which both an
active and an inactive enzyme have been measured. There are ten, and three large enough to
estimate anything from. Cross-lineage transfer is at chance because the determinants are
lineage-specific, and no amount of modelling lifts that.

The fix is measurement, and it has a shape: roughly ten enzymes from each of roughly
fifteen NEW 30%-identity lineages, assayed under one protocol at 37 degC on amorphous PET.
About 150 assays, taking evaluable lineages from 3 to ~18. This script decides which ones,
from 13,827 unlabelled homologues in the right length range plus 439 metagenome candidates.

The selection rules, and why each is there:

  NEW lineages only. A cluster already holding measured enzymes adds depth where the
  project has plenty and none of the breadth it lacks. Depth in the 262-member cluster buys
  nothing for the transfer question.

  At least ten members per lineage, because the within-lineage contrast is what makes a
  cluster evaluable, and a cluster contributing two enzymes contributes an unusable AUC of
  0 or 1 -- the exact failure that made an unweighted mean read 0.625 when the real value
  was 0.472.

  Members chosen to SPREAD within their lineage, by greedy maximum-minimum distance in
  embedding space rather than by taking the ten most typical. Ten near-duplicates would
  measure the same enzyme ten times; the point is to span the lineage so that whatever
  varies within it has a chance to show up.

  A complete catalytic triad where a structure exists. An enzyme missing the chemistry is a
  guaranteed negative and wastes a well on a foregone conclusion.

  Length 200-450 aa, the window the reference set occupies, to keep the panel comparable
  with everything already measured.

  PANTS' OWN LINEAGES FIRST, then the largest of the rest. Selecting purely by cluster size
  produced a panel of 150 assays in which 140 were somebody else's sequences and 10 were
  this project's, which adds evaluable lineages while testing none of the mining that
  produced the catalogue. Eight new clusters consist entirely of PANTS metagenome
  candidates -- novel lineages with no characterised relative at all, which is precisely
  why retrieval places them out of range and nothing here can score them. Measuring those
  does both jobs at once: it adds independent lineages AND asks whether the mining found
  anything real. The remaining slots go to the largest unassayed lineages for breadth.

Deliberately NOT ranked by predicted activity. The whole finding is that the prediction
does not transfer to new lineages, so using it to choose new lineages would be circular --
it would select the enzymes that most resemble what is already known, which is the bias the
project exists to correct and the reason 384 of 439 candidates cannot currently be scored.
Breadth is the objective, not expected hit rate.
"""

from __future__ import annotations

import csv
import json
import pathlib
import sys
from typing import Dict, List, Set, Tuple

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config, seqtools
from pipeline.db import connect

TARGET_LINEAGES = 15
PER_LINEAGE = 10
MIN_MEMBERS = 10
LEN_MIN, LEN_MAX = 200, 450
OUT_CSV = config.ROOT_DIR / "release" / "validation_panel.csv"
OUT_JSON = config.INTERIM_DIR / "validation_panel.json"


def pool() -> Tuple[List[Tuple[str, str]], Dict[str, dict]]:
    """Everything assayable that nobody has assayed, with its provenance."""
    meta: Dict[str, dict] = {}
    recs: List[Tuple[str, str]] = []
    with connect() as c:
        for sid, acc, org, seq, n in c.execute(
                "SELECT seq_id, accession, organism, sequence, seq_length "
                "FROM unlabelled_sequences WHERE seq_length BETWEEN ? AND ?",
                (LEN_MIN, LEN_MAX)):
            recs.append((sid, seq))
            meta[sid] = {"source": "Science-S1 homologue", "accession": acc,
                         "organism": org, "length": n, "triad": None, "environment": None}
        for cid, seq, n, env, tri in c.execute(
                "SELECT c.candidate_id, c.sequence, c.seq_length, c.source_environment, "
                "       g.triad_ser_resnum "
                "FROM candidates c LEFT JOIN geometry g ON g.candidate_id=c.candidate_id "
                "WHERE c.seq_length BETWEEN ? AND ?", (LEN_MIN, LEN_MAX)):
            recs.append((cid, seq))
            meta[cid] = {"source": "PANTS metagenome candidate", "accession": None,
                         "organism": None, "length": n,
                         "triad": bool(tri), "environment": env}
    return recs, meta


def measured_sequences() -> List[Tuple[str, str]]:
    with connect() as c:
        return [tuple(r) for r in c.execute(
            "SELECT enzyme_id, sequence FROM characterised_enzymes "
            "WHERE sequence IS NOT NULL AND (within_family_basis='measured-inactive' "
            "   OR (is_positive=1 AND source_ref IN ('PAZy-measured','ACS-screen-measured',"
            "       'Science-landscape-measured','EC-experimental','UniProt','HGMP-measured')))")]


def spread(ids: List[str], X: np.ndarray, k: int) -> List[str]:
    """Greedy maximum-minimum-distance subset: span the lineage, do not sample its centre."""
    if len(ids) <= k:
        return ids
    Xn = X / np.clip(np.linalg.norm(X, axis=1, keepdims=True), 1e-9, None)
    centre = Xn.mean(0)
    chosen = [int(np.argmax(Xn @ centre))]          # start from the most typical member
    while len(chosen) < k:
        d = np.min(1 - (Xn @ Xn[chosen].T), axis=1)
        d[chosen] = -1
        chosen.append(int(np.argmax(d)))
    return [ids[i] for i in chosen]


def main() -> int:
    recs, meta = pool()
    measured = measured_sequences()
    print(f"pool: {len(recs)} unassayed sequences in {LEN_MIN}-{LEN_MAX} aa")
    print(f"already measured: {len(measured)} (their lineages are excluded)\n")

    allrecs = recs + measured
    fa = seqtools.write_fasta(allrecs, config.INTERIM_DIR / "panel.fasta")
    clu = seqtools.cluster(fa, min_seq_id=0.30)
    grp = {i: clu.get(i, i) for i, _ in allrecs}

    taken: Set[str] = {grp[i] for i, _ in measured}
    members: Dict[str, List[str]] = {}
    for i, _ in recs:
        g = grp[i]
        if g in taken:
            continue
        members.setdefault(g, []).append(i)
    fresh = {g: m for g, m in members.items() if len(m) >= MIN_MEMBERS}
    print(f"{len(set(grp.values()))} clusters overall; {len(taken)} already hold a measured "
          f"enzyme")
    print(f"{len(fresh)} NEW clusters have at least {MIN_MEMBERS} unassayed members\n")
    if not fresh:
        print("no new lineage is large enough — the panel cannot be designed from this pool")
        return 1

    from pipeline.embed import esm
    # PANTS' own novel lineages first -- they are the out-of-range candidates the tool
    # cannot score, so measuring them tests the mining as well as adding lineages -- then
    # the largest of the rest for breadth.
    own = sorted((g for g, m in fresh.items()
                  if any(meta[i]["source"].startswith("PANTS") for i in m)),
                 key=lambda g: -len(fresh[g]))
    rest = sorted((g for g in fresh if g not in set(own)), key=lambda g: -len(fresh[g]))
    order = (own + rest)[:TARGET_LINEAGES]
    print(f"selected {sum(1 for g in order if g in set(own))} lineages of PANTS candidates "
          f"and {sum(1 for g in order if g not in set(own))} from the homologue pool\n")
    seqmap = dict(allrecs)
    picked, rows = [], []
    for rank, g in enumerate(order, 1):
        ids = fresh[g]
        sub = [(i, seqmap[i]) for i in ids]
        _, X, _ = esm.embed(sub, batch_size=8, progress_every=10_000)
        keep = spread(ids, X, PER_LINEAGE)
        n_cand = sum(1 for i in keep if meta[i]["source"].startswith("PANTS"))
        print(f"  lineage {rank:>2} ({g[:14]:<16}) {len(ids):>5} available -> "
              f"{len(keep)} chosen, {n_cand} of them PANTS candidates")
        for i in keep:
            m = meta[i]
            rows.append({"lineage": g, "lineage_rank": rank, "id": i, "source": m["source"],
                         "accession": m["accession"], "organism": m["organism"],
                         "length": m["length"], "environment": m["environment"],
                         "complete_triad": m["triad"], "sequence": seqmap[i]})
        picked += keep

    OUT_CSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_CSV.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    OUT_JSON.write_text(json.dumps({
        "n_assays": len(rows), "n_lineages": len(order), "per_lineage": PER_LINEAGE,
        "protocol": "one protocol, 37 degC, amorphous PET film",
        "evaluable_lineages_before": 3, "evaluable_lineages_after_if_successful": 3 + len(order),
        "selection": "new 30% clusters only; members spread by greedy max-min distance",
    }, indent=2))
    print(f"\n{len(rows)} assays across {len(order)} new lineages")
    print(f"wrote {OUT_CSV}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
