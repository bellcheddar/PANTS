"""SQLite persistence: shared by the pipeline (writer) and the Flask app (reader).

Connection handling is ported from AlphaFraud's db.py, including two things it learned in
production and that are easy to get wrong:

  1. WAL is a PERSISTENT property of the database file. Set it once in init_schema(),
     never on every connection: `PRAGMA journal_mode=WAL` needs a brief write lock, and
     when a deploy restarts the web service (which also has the DB open) that collides
     with pipeline writes and surfaces as "attempt to write a readonly database".
  2. Wrap writes in retry_write(). A transient OperationalError from a locked window
     should not crash a multi-hour batch stage.
"""

from __future__ import annotations

import sqlite3
import time
from contextlib import contextmanager
from datetime import datetime, timezone
from typing import Any, Callable, Iterator, Optional

from .. import config
from .schema import COLUMN_MIGRATIONS, SCHEMA, SCHEMA_VERSION

__all__ = [
    "connect", "init_schema", "retry_write", "now",
    "start_run", "finish_run", "fail_run",
    "SCHEMA_VERSION",
]


def now() -> str:
    """UTC timestamp, ISO 8601, second resolution. One format everywhere in the DB."""
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


@contextmanager
def connect() -> Iterator[sqlite3.Connection]:
    conn = sqlite3.connect(config.DB_PATH, timeout=30)
    conn.row_factory = sqlite3.Row
    # Per-connection and harmless: wait for a writer rather than erroring immediately.
    conn.execute("PRAGMA busy_timeout=30000")
    # Foreign keys are OFF by default in SQLite and must be enabled per connection.
    conn.execute("PRAGMA foreign_keys=ON")
    try:
        yield conn
        conn.commit()
    finally:
        conn.close()


def init_schema() -> None:
    """Create the schema if missing. Safe to call repeatedly (every statement is IF NOT EXISTS)."""
    with connect() as conn:
        conn.execute("PRAGMA journal_mode=WAL")   # persisted in the DB header; set once
        conn.executescript(SCHEMA)
        _apply_column_migrations(conn)
        conn.execute(
            "INSERT INTO app_state(key, value) VALUES ('schema_version', ?) "
            "ON CONFLICT(key) DO UPDATE SET value=excluded.value",
            (str(SCHEMA_VERSION),),
        )


def _apply_column_migrations(conn) -> None:
    """Add columns that the CREATE TABLE IF NOT EXISTS statements cannot reach.

    See COLUMN_MIGRATIONS in schema.py: on an existing database the CREATE is skipped
    entirely, so a newly declared column is simply absent until it is ALTERed in.
    """
    for table, column, decl in COLUMN_MIGRATIONS:
        cols = {r[1] for r in conn.execute(f"PRAGMA table_info({table})")}
        if not cols:
            continue                      # table itself does not exist yet
        if column not in cols:
            conn.execute(f"ALTER TABLE {table} ADD COLUMN {column} {decl}")


def retry_write(action: Callable[[], Any], attempts: int = 6, base_delay: float = 0.4) -> Any:
    """Retry a DB write through a transient OperationalError (a locked or readonly window,
    e.g. a concurrent web-service restart) instead of letting a long batch stage crash."""
    for i in range(attempts):
        try:
            return action()
        except sqlite3.OperationalError:
            if i == attempts - 1:
                raise
            time.sleep(base_delay * (i + 1))


# --------------------------------------------------------------------------------------
# Runs: every pipeline stage opens one, and every manifest hangs off it
# --------------------------------------------------------------------------------------
def start_run(stage: str, label: str, params_json: Optional[str] = None) -> int:
    def _do() -> int:
        with connect() as conn:
            cur = conn.execute(
                "INSERT INTO runs(stage, label, started_at, status, params_json) "
                "VALUES (?,?,?, 'running', ?)",
                (stage, label, now(), params_json),
            )
            return int(cur.lastrowid)
    return retry_write(_do)


def finish_run(run_id: int, n_input: int = 0, n_output: int = 0, n_discarded: int = 0) -> None:
    """Close a run as done. n_discarded matters: spec section 5.1 requires the recall stage
    to report how many candidates were dropped at the triad filter."""
    def _do() -> None:
        with connect() as conn:
            conn.execute(
                "UPDATE runs SET finished_at=?, n_input=?, n_output=?, n_discarded=?, "
                "status='done' WHERE id=?",
                (now(), n_input, n_output, n_discarded, run_id),
            )
    retry_write(_do)


def fail_run(run_id: int) -> None:
    def _do() -> None:
        with connect() as conn:
            conn.execute(
                "UPDATE runs SET finished_at=?, status='error' WHERE id=?", (now(), run_id)
            )
    retry_write(_do)
