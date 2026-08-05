"""Assemble the deposit bundle, and derive every published number from the artefacts.

The verification half is not ceremony. On the sibling TopPDBLX deposit, six figures in the
Zenodo description were wrong because they had been copied from README prose describing a
checkpoint that was measured but never shipped, and one had already been pasted into a
description where a DOI would have made it permanent. The lesson recorded there:
**verify published numbers against artefacts, never against your own earlier prose.**

So this script exports the tables, then recomputes every headline figure from the exported
files themselves and writes them to release/STATS.json. The Zenodo description is built
from that JSON, so a number can only appear if it survived a round trip through the data.
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

MEASURED_TIERS = ("EC-experimental", "UniProt", "HGMP-measured", "PAZy-measured",
                  "Tournier et al. 2020, Nature", "Son et al. 2019, ACS Catal.",
                  "Lu et al. 2022, Nature (MutCompute)", "Austin et al. 2018, PNAS")

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
