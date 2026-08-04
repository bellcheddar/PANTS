"""Phase 0 tests: schema, run lifecycle, manifest provenance.

Every test runs against a throwaway database via the PANTS_DB env var, so the real
pants.db is never touched.
"""

from __future__ import annotations

import importlib
import json

import pytest


@pytest.fixture()
def pants(tmp_path, monkeypatch):
    """A freshly initialised PANTS database in tmp_path, with config/db reloaded to see it."""
    monkeypatch.setenv("PANTS_DB", str(tmp_path / "test.db"))
    from pipeline import config as _config
    importlib.reload(_config)
    from pipeline import db as _db
    importlib.reload(_db)
    from pipeline.db import manifest as _manifest
    importlib.reload(_manifest)

    monkeypatch.setattr(_config, "MANIFEST_DIR", tmp_path / "manifests")
    monkeypatch.setattr(_manifest.config, "MANIFEST_DIR", tmp_path / "manifests")
    _db.init_schema()
    return _db, _manifest, _config


def test_schema_creates_every_table(pants):
    db, _, _ = pants
    expected = {
        "runs", "manifests", "data_sources", "candidates", "scores", "structures",
        "geometry", "characterised_enzymes", "activity_measurements", "eval_splits",
        "training_runs", "app_state", "visits",
    }
    with db.connect() as conn:
        got = {r[0] for r in conn.execute(
            "SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'"
        )}
    assert expected <= got


def test_init_schema_is_idempotent(pants):
    db, _, _ = pants
    db.init_schema()
    db.init_schema()
    with db.connect() as conn:
        assert conn.execute(
            "SELECT value FROM app_state WHERE key='schema_version'"
        ).fetchone()["value"] == str(db.SCHEMA_VERSION)


def test_wal_is_enabled(pants):
    db, _, _ = pants
    with db.connect() as conn:
        assert conn.execute("PRAGMA journal_mode").fetchone()[0].lower() == "wal"


def test_run_lifecycle_records_discard_count(pants):
    """spec section 5.1 requires the recall stage to report how many candidates it drops."""
    db, _, _ = pants
    run_id = db.start_run("recall", "unit-test")
    db.finish_run(run_id, n_input=1000, n_output=40, n_discarded=960)
    with db.connect() as conn:
        row = conn.execute("SELECT * FROM runs WHERE id=?", (run_id,)).fetchone()
    assert row["status"] == "done"
    assert (row["n_input"], row["n_output"], row["n_discarded"]) == (1000, 40, 960)


def test_manifest_records_hashes_and_versions(pants, tmp_path):
    db, manifest, config = pants
    src = tmp_path / "in.fasta"
    src.write_text(">a\nMKV\n")
    out = tmp_path / "hits.tsv"

    with manifest.stage_manifest("recall", label="unit", inputs=[src]) as m:
        out.write_text("a\t1e-40\n")
        m.add_output(out)
        m.counts(n_input=1, n_output=1, n_discarded=0)

    with db.connect() as conn:
        row = conn.execute("SELECT * FROM manifests ORDER BY id DESC LIMIT 1").fetchone()

    assert row["stage"] == "recall"
    assert row["schema_version"] == db.SCHEMA_VERSION
    assert json.loads(row["input_hashes_json"])[str(src)] != "missing"
    assert json.loads(row["output_hashes_json"])[str(out)] != "missing"
    # mmseqs2 and hmmer are both installed, so both versions must be captured.
    tools = json.loads(row["tool_versions_json"])
    assert tools["mmseqs2"]
    assert tools["hmmer"] and "HMMER" in tools["hmmer"], \
        f"hmmer version should carry the version string, got {tools['hmmer']!r}"


def test_failed_stage_still_leaves_provenance(pants):
    """A stage that raises must mark the run 'error' AND write its manifest, so a failure
    leaves a trail rather than nothing."""
    db, manifest, _ = pants
    with pytest.raises(RuntimeError):
        with manifest.stage_manifest("embed", label="deliberate-failure"):
            raise RuntimeError("boom")

    with db.connect() as conn:
        run = conn.execute("SELECT * FROM runs ORDER BY id DESC LIMIT 1").fetchone()
        n = conn.execute("SELECT COUNT(*) FROM manifests WHERE stage='embed'").fetchone()[0]
    assert run["status"] == "error"
    assert n == 1


def test_sha256_flags_partial_hashes(pants, tmp_path, monkeypatch):
    """Large files get a head+tail digest; it must be labelled so it is never mistaken for
    a full one."""
    _, manifest, _ = pants
    monkeypatch.setattr(manifest, "FULL_HASH_LIMIT_BYTES", 16)
    big = tmp_path / "big.fasta"
    big.write_bytes(b"M" * 4096)
    assert manifest.sha256_file(big).startswith("partial:")

    small = tmp_path / "small.fasta"
    small.write_bytes(b"M" * 8)
    assert not manifest.sha256_file(small).startswith("partial:")


def test_missing_file_hashes_to_missing(pants, tmp_path):
    _, manifest, _ = pants
    assert manifest.sha256_file(tmp_path / "nope.fasta") == "missing"
