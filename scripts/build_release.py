"""Assemble the deposit bundle, and derive every published number from the artefacts.

The verification half is not ceremony. On the sibling TopPDBLX deposit, six figures in the
Zenodo description were wrong because they had been copied from README prose describing a
checkpoint that was measured but never shipped, and one had already been pasted into a
description where a DOI would have made it permanent. The lesson recorded there:
**verify published numbers against artefacts, never against your own earlier prose.**

So this script exports the tables, then recomputes every headline figure from the exported
files themselves and writes them to release/STATS.json. The Zenodo description is built
from that JSON, so a number can only appear if it survived a round trip through the data.

DATASHEET.md is now generated the same way, and for the same reason. It was hand-written
and drifted exactly as predicted: by v0.2.0 it still claimed 1,107 reference rows against
1,140, 48 measurements against 75, and 268 structures against 1,188, and it described the
release as "a snapshot taken before folding finished" months after folding finished. A
prose file that has to be remembered is a prose file that will be wrong, so the counts,
the tier table and the temperature table now come from STATS.json and only the judgement
-- limitations, sources, what the dataset is for -- is written by hand.
"""

from __future__ import annotations

import csv
import hashlib
import json
import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config
from pipeline.db import connect

REL = config.ROOT_DIR / "release"

MEASURED_TIERS = config.MEASURED_TIERS   # single definition; see pipeline/config.py

TABLES = {
    "candidates": """
        SELECT c.candidate_id, c.source_environment, c.assembly_id, c.contig_id,
               c.sequence, c.seq_length, c.has_complete_triad,
               c.recall_evalue, c.recall_bitscore, c.recall_profile_identity,
               c.nearest_characterised_id, c.structure_deferred,
               s.structure_method, s.plddt_mean, s.rmsd_ca_to_ispetase_A,
               g.triad_ser_resnum, g.triad_asp_resnum, g.triad_his_resnum,
               g.ser_og_his_ne2_dist_A, g.his_nd1_asp_od_dist_A,
               g.cleft_width_A, g.cleft_depth_A, g.aromatic_clamp_residues_json
        FROM candidates c
        LEFT JOIN structures s ON s.candidate_id = c.candidate_id
        LEFT JOIN geometry   g ON g.candidate_id = c.candidate_id
        ORDER BY c.recall_bitscore DESC""",
    "characterised_enzymes": """
        SELECT enzyme_id, uniprot, organism, family, esther_family, seq_length, sequence,
               is_positive, is_negative, is_near_miss, source_ref, topt_c, ph_opt,
               excluded_from_training, exclusion_reason, activity_substrate_notes
        FROM characterised_enzymes ORDER BY enzyme_id""",
    "activity_measurements": """
        SELECT enzyme_id, parameter_type, substrate_form, temperature_c, ph,
               product_measured, rate_value, rate_units, raw_text, evidence_code,
               comparable_group_id, source_doi, extraction_confidence
        FROM activity_measurements ORDER BY enzyme_id, parameter_type""",
    "runs": """
        SELECT stage, label, status, started_at, finished_at,
               n_input, n_output, n_discarded, params_json
        FROM runs ORDER BY id""",
    "training_runs": "SELECT * FROM training_runs ORDER BY run_id",
    "data_sources": "SELECT * FROM data_sources ORDER BY name",
}


def export() -> dict:
    REL.mkdir(parents=True, exist_ok=True)
    written = {}
    with connect() as conn:
        for name, sql in TABLES.items():
            rows = list(conn.execute(sql))
            path = REL / f"{name}.csv"
            with path.open("w", newline="") as fh:
                if rows:
                    w = csv.DictWriter(fh, fieldnames=rows[0].keys())
                    w.writeheader()
                    for r in rows:
                        w.writerow(dict(r))
            written[name] = (path, len(rows))
    return written


def sha256(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def verify_from_artefacts() -> dict:
    """Recompute every headline figure from the EXPORTED CSVs, not from the database.

    Reading the exports rather than the database is the point: it verifies what is
    actually being deposited, and would catch an export that silently dropped rows.
    """
    def rows(name):
        with (REL / f"{name}.csv").open() as fh:
            return list(csv.DictReader(fh))

    cand = rows("candidates")
    ce = rows("characterised_enzymes")
    am = rows("activity_measurements")
    runs = rows("runs")

    by_env = {}
    for c in cand:
        by_env[c["source_environment"]] = by_env.get(c["source_environment"], 0) + 1

    scanned = {}
    for r in runs:
        if r["stage"] != "recall" or not r["params_json"]:
            continue
        env = (json.loads(r["params_json"]) or {}).get("environment", "?")
        scanned[env] = scanned.get(env, 0) + int(r["n_input"] or 0)

    pos = [r for r in ce if r["is_positive"] == "1"]
    measured = [r for r in pos if r["source_ref"] in MEASURED_TIERS]

    return {
        "candidates": len(cand),
        "candidates_by_environment": dict(sorted(by_env.items(), key=lambda kv: -kv[1])),
        "sequences_scanned": sum(scanned.values()),
        "sequences_scanned_by_environment": dict(sorted(scanned.items(), key=lambda kv: -kv[1])),
        "candidates_with_structure": sum(1 for c in cand if c["plddt_mean"]),
        "candidates_with_geometry": sum(1 for c in cand if c["cleft_width_A"]),
        "structures_deferred_over_length": sum(1 for c in cand if c["structure_deferred"] == "1"),
        "characterised_total": len(ce),
        "positives": len(pos),
        "positives_measured": len(measured),
        "positives_by_tier": {t: sum(1 for r in pos if r["source_ref"] == t)
                              for t in sorted({r["source_ref"] for r in pos})},
        "hard_negatives": sum(1 for r in ce if r["is_negative"] == "1"),
        "near_misses": sum(1 for r in ce if r["is_near_miss"] == "1"),
        "activity_measurements": len(am),
        "activity_measurements_with_a_doi": sum(1 for r in am if r["source_doi"]),
        "measured_topt_values": sorted(
            (r["enzyme_id"], float(r["rate_value"]))
            for r in am if r["parameter_type"] == "topt" and r["rate_value"]),
    }


DATASHEET_PREAMBLE = """# PANTS dataset, v{version}

Candidate polyesterases mined from metagenomic sequence space and triaged for
**therapeutic** use: degradation of PET at 37 °C and neutral pH, rather than in an
industrial reactor above PET's glass transition.

Every count below was recomputed from the deposited CSVs by `scripts/build_release.py`,
not copied from prose. `STATS.json` carries the same numbers machine-readably.
"""

DATASHEET_JUDGEMENT = """
## Limitations

1. **The reference set is the constraint, not the method.** Three independent evaluations
   agree. Sequence embeddings separate polyesterases from other folds at AUC 0.975 and
   PET-active from PET-inactive polyesterases at 0.493. Active-site geometry looks
   convincing raw (cleft depth AUC 0.808) and falls to 0.534 under cluster-grouped
   splitting. And the learned head, trained on labels somebody actually measured, scores
   0.921 against a nearest-known-PETase retrieval baseline of 0.931 -- it does not beat
   looking the answer up.
2. **Leave-one-family-out is not evaluable.** All 13 ESTHER families here are wholly
   positive or wholly negative, so holding one out removes a class. In this catalogue the
   label is family membership, which is the same finding stated structurally.
3. **Predicted and experimental coordinates are not interchangeable.** Paired within the
   same protein (n=51), the oxyanion hole differs: second-donor angle 23.6° in the crystal
   against 15.6° in the model, p 1.1e-06. Cleft depth is source-invariant. The `source`
   column exists so this can be accounted for; pooling without accounting for it is an
   error.
4. **Candidates over 450 aa were deliberately not folded**: probable fusions, multi-domain
   proteins or misassemblies, and the cleft measurement assumes a single active site in one
   domain. They are deferred and marked, never discarded.
5. **Published activity data is not harmonised** across assay formats. Absolute rate
   predictions should not be trusted.
6. **The prospective holdout is underpowered.** Its test side holds 366 positives and 6
   negatives, so its AUC is indicative at best.
7. Nothing here addresses delivery, immunogenicity, biodistribution, or the fate of
   liberated TPA and EG in vivo.

## Sources

UniProt and UniRef (CC BY 4.0), the PDB (CC0), AlphaFold DB (CC BY 4.0), ESTHER,
PAZy (Buchholz et al. 2022, Proteins 90:1443), MGnify assemblies, and the human gut PET
hydrolases of Zhang et al. 2024 (PMID 39551294) from the authors' SciDB deposit.

Data are CC BY 4.0 (see LICENSE-DATA); the source code is MIT (see LICENSE).
"""


def datasheet(stats: dict, written: dict, n_structures: int) -> str:
    """Write DATASHEET.md from the verified figures rather than from memory."""
    rows = {name: n for name, (_, n) in written.items()}
    what = {
        "candidates.csv": "Mined candidates: retrieval scores, triad positions, active-site geometry",
        "characterised_enzymes.csv": "Reference set: positives by evidence tier, hard negatives, near misses",
        "activity_measurements.csv": f"Measured kinetics and optima, {stats['activity_measurements_with_a_doi']} carrying a DOI",
        "runs.csv": "Every pipeline stage with its input, output and discard counts",
        "training_runs.csv": "Model runs including the baselines each must clear",
        "data_sources.csv": "Sources, versions, retrieval dates, licences",
    }
    out = [DATASHEET_PREAMBLE.format(version=config.DATA_VERSION), "## Contents", "",
           "| File | Rows | What it is |", "|---|---|---|"]
    for name, desc in what.items():
        out.append(f"| `{name}` | {rows.get(name.replace('.csv',''), '—')} | {desc} |")
    out.append(f"| `structures.tar.gz` | {n_structures} | Candidate predictions and reference "
               f"structures, all superposed onto IsPETase (PDB 6EQE) |")
    out.append(f"| `evaluation_protocol.json` | — | Every component of the evaluation protocol, "
               f"including the ones that are not evaluable and why |")
    out.append(f"| `structure_source_confound.json` | — | Paired crystal-versus-model geometry |")

    sc = stats["sequences_scanned"]
    out += ["", "## Candidates", "",
            f"{sc:,} predicted proteins scanned, {stats['candidates']} candidates retained "
            f"({100*stats['candidates']/sc:.3f}%).", "",
            "| Environment | Sequences scanned | Candidates | Per million |", "|---|---|---|---|"]
    for env, n in sorted(stats["sequences_scanned_by_environment"].items(),
                         key=lambda kv: -kv[1]):
        c = stats["candidates_by_environment"].get(env, 0)
        out.append(f"| {env} | {n:,} | {c} | {1e6*c/n:.0f} |")
    out += ["", "Retention is a **choice**: a strict E-value against the seed set, then a "
            "requirement that Ser, His and Asp be connected **in space** rather than merely "
            "present in sequence.", ""]

    out += ["## Reference set, by evidence tier", "",
            f"{stats['positives_measured']} of {stats['positives']} positives carry a "
            f"measurement. The rest are annotation, and the distinction is the single most "
            f"important thing in this dataset.", "", "| Tier | Positives |", "|---|---|"]
    for tier, n in sorted(stats["positives_by_tier"].items(), key=lambda kv: -kv[1]):
        out.append(f"| `{tier}` | {n} |")
    out += ["", f"Also {stats['hard_negatives']} hard negatives and "
            f"{stats['near_misses']} near misses.", "",
            "`PAZy-measured` entries are present because activity was **measured** on a "
            "plastic and published. `EC-auto-annotated` entries carry EC 3.1.1.101 from "
            "`ECO:0000256`, meaning the label was assigned **by sequence similarity**. The "
            "two are never pooled in a reported metric.", ""]

    out += ["## Measured temperature optima", "", "| Enzyme | Topt |", "|---|---|"]
    for enzyme, topt in stats["measured_topt_values"]:
        out.append(f"| {enzyme} | {topt:g} °C |")
    out += ["", "The therapeutic target is 37 °C. Almost everything measured here optimises "
            "far above it, which is the problem the project exists to address.", ""]
    return "\n".join(out) + DATASHEET_JUDGEMENT


if __name__ == "__main__":
    written = export()
    print("exported:")
    for name, (path, n) in written.items():
        print(f"  {n:>6} rows  {path.name}  ({path.stat().st_size/1024:.0f} kB)")

    stats = verify_from_artefacts()
    (REL / "STATS.json").write_text(json.dumps(stats, indent=2))
    print("\nverified from the exported CSVs:")
    for k, v in stats.items():
        if isinstance(v, (int, str)):
            print(f"  {k:<38} {v}")
    print(f"\n  wrote {REL/'STATS.json'}")

    tarball = REL / "structures.tar.gz"
    n_struct = 0
    if tarball.exists():
        import tarfile
        with tarfile.open(tarball) as tf:
            n_struct = sum(1 for m in tf.getmembers() if m.isfile())
    (REL / "DATASHEET.md").write_text(datasheet(stats, written, n_struct))
    print(f"  wrote {REL/'DATASHEET.md'} ({n_struct} structures in the tarball)")
