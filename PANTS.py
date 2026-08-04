#!/usr/bin/env python3
"""PANTS: PETase ANnotation and Triage System.

Command line entry point for the offline pipeline and the local dev server, mirroring
AlphaFraud.py's init/run/status/serve shape.

  ./PANTS.py init      create dirs, DB schema, check the data symlinks
  ./PANTS.py status    database summary: candidates, scores, structures, last runs
  ./PANTS.py serve     local dev web server (production uses gunicorn, see deploy/)

Pipeline subcommands (recall, embed, train, structure) are added as each stage lands.
"""

from __future__ import annotations

import argparse
import sys

from pipeline import config
from pipeline import db


def cmd_init(_args) -> int:
    config.ensure_dirs()
    db.init_schema()
    print(f"database:  {config.DB_PATH}")
    print(f"schema:    v{db.SCHEMA_VERSION}")
    print(f"manifests: {config.MANIFEST_DIR}")

    if config.data_dirs_are_external():
        print(f"data/raw:  -> {config.RAW_DIR.resolve()}  (outside iCloud, correct)")
    else:
        print("WARNING: data/raw and data/interim are NOT symlinks out of the repo.")
        print("         Documents is iCloud-synced and macOS 'Optimize Mac Storage' will")
        print("         evict large files mid-run. Fix with:")
        print("           mkdir -p ~/PANTSData/{raw,interim}")
        print("           rmdir data/raw data/interim")
        print("           ln -s ~/PANTSData/raw data/raw")
        print("           ln -s ~/PANTSData/interim data/interim")
    return 0


def cmd_status(_args) -> int:
    with db.connect() as conn:
        def count(table: str, where: str = "") -> int:
            try:
                return int(conn.execute(f"SELECT COUNT(*) FROM {table} {where}").fetchone()[0])
            except Exception:
                return -1

        print("PANTS database:", config.DB_PATH)
        print()
        print(f"  candidates            {count('candidates'):>8}")
        print(f"    complete triad      {count('candidates', 'WHERE has_complete_triad=1'):>8}")
        print(f"  scored                {count('scores'):>8}")
        print(f"  structures            {count('structures'):>8}")
        print(f"  geometry rows         {count('geometry'):>8}")
        print(f"  characterised         {count('characterised_enzymes'):>8}")
        print(f"    positives           {count('characterised_enzymes', 'WHERE is_positive=1'):>8}")
        print(f"    hard negatives      {count('characterised_enzymes', 'WHERE is_negative=1'):>8}")
        print(f"    near misses         {count('characterised_enzymes', 'WHERE is_near_miss=1'):>8}")
        print(f"  activity measurements {count('activity_measurements'):>8}")
        print(f"  training runs         {count('training_runs'):>8}")

        rows = conn.execute(
            "SELECT stage, label, status, started_at, n_input, n_output, n_discarded "
            "FROM runs ORDER BY id DESC LIMIT 8"
        ).fetchall()
        if rows:
            print("\n  recent runs:")
            for r in rows:
                print(f"    {r['started_at']}  {r['stage']:<10} {r['status']:<8} "
                      f"in={r['n_input']} out={r['n_output']} discarded={r['n_discarded']}"
                      f"  ({r['label']})")
        else:
            print("\n  no runs yet")
    return 0


def cmd_curate_seeds(args) -> int:
    """Fetch the characterised wild types from UniProt and derive the engineered variants."""
    from pipeline.recall.load import load_reference_set

    report = load_reference_set(label=args.label)

    print("wild types (fetched from UniProt):")
    for enzyme_id, length, organism in report["wild_types"]:
        print(f"  {enzyme_id:<22} {length:>4} aa   {organism}")

    print("\nvariants (derived from parent, mutations validated):")
    for enzyme_id, status, offset, length in report["variants"]:
        detail = f"{length} aa, offset {offset:+d}" if status == "derived" else status
        print(f"  {enzyme_id:<22} {detail}")

    if report["problems"]:
        print("\nPROBLEMS:")
        for p in report["problems"]:
            print(f"  {p}")
        return 1
    return 0


def cmd_serve(args) -> int:
    from app import create_app
    create_app().run(host=args.host, port=args.port, debug=args.debug)
    return 0


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="PANTS", description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    sub = parser.add_subparsers(dest="command", required=True)

    sub.add_parser("init", help="create dirs, DB schema, check data symlinks"
                   ).set_defaults(func=cmd_init)
    sub.add_parser("status", help="database summary"
                   ).set_defaults(func=cmd_status)

    p_cs = sub.add_parser("curate-seeds",
                          help="fetch characterised wild types from UniProt, derive variants")
    p_cs.add_argument("--label", default="v1", help="run label recorded in the manifest")
    p_cs.set_defaults(func=cmd_curate_seeds)

    p_serve = sub.add_parser("serve", help="local dev web server")
    p_serve.add_argument("--host", default="127.0.0.1")
    p_serve.add_argument("--port", type=int, default=8005)
    p_serve.add_argument("--debug", action="store_true")
    p_serve.set_defaults(func=cmd_serve)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
