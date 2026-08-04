"""Manifest writing: the provenance record every pipeline stage must leave behind.

Spec section 12: "Every pipeline stage writes a manifest with input hashes, tool versions
and model version. The app displays the current data version in the footer."

A manifest is written twice, deliberately:
  - as a row in the `manifests` table, so the Methods tab can render it without touching
    the filesystem, and
  - as a JSON file under manifests/, so provenance survives the database being rebuilt.

Usage, from any stage:

    with stage_manifest("recall", label="mgnify-MGYS00001234",
                        inputs=[fasta_path], tool_versions=tool_versions()) as m:
        ...do the work...
        m.add_output(hits_path)
        m.counts(n_input=3_000_000, n_output=4_812, n_discarded=41_009)
"""

from __future__ import annotations

import hashlib
import json
import subprocess
import time
from contextlib import contextmanager
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, Iterable, Iterator, Optional

from .. import config
from . import connect, fail_run, finish_run, now, retry_write, start_run
from .schema import SCHEMA_VERSION

# Hashing a multi-GB metagenome FASTA in full costs minutes for no extra provenance value,
# so files above this size are hashed over head+tail+size instead. The manifest records
# which method was used so the number is never mistaken for a full digest.
FULL_HASH_LIMIT_BYTES = 256 * 1024 * 1024
_CHUNK = 1024 * 1024


def sha256_file(path: str | Path) -> str:
    """SHA256 of a file. Large files get a head+tail+size digest, prefixed 'partial:'."""
    p = Path(path)
    if not p.exists():
        return "missing"
    size = p.stat().st_size
    h = hashlib.sha256()
    if size <= FULL_HASH_LIMIT_BYTES:
        with p.open("rb") as fh:
            for chunk in iter(lambda: fh.read(_CHUNK), b""):
                h.update(chunk)
        return h.hexdigest()
    # Clamp the tail read: seeking -_CHUNK from the end runs off the front of the file
    # whenever size < _CHUNK, which raises OSError(EINVAL) rather than reading less.
    tail_len = min(_CHUNK, size)
    with p.open("rb") as fh:
        head = fh.read(_CHUNK)
        fh.seek(-tail_len, 2)
        tail = fh.read(tail_len)
    h.update(head)
    h.update(tail)
    h.update(str(size).encode())
    return "partial:" + h.hexdigest()


def _cmd_version(args: list[str], match: Optional[str] = None) -> Optional[str]:
    """First line of a tool's version output, or the first line containing `match`.

    HMMER needs the `match` form: `hmmbuild -h` opens with a usage banner and only puts
    the version on a later line ("# HMMER 3.4 (Aug 2023)"), so taking line 0 records the
    banner as the version and the manifest silently loses the number that matters.
    """
    try:
        out = subprocess.run(args, capture_output=True, text=True, timeout=15)
        lines = [ln.strip() for ln in (out.stdout + out.stderr).splitlines() if ln.strip()]
        if not lines:
            return None
        if match:
            hit = next((ln for ln in lines if match in ln), None)
            return hit.lstrip("# ").strip() if hit else None
        return lines[0]
    except (OSError, subprocess.SubprocessError):
        return None


def git_commit() -> Optional[str]:
    try:
        out = subprocess.run(
            ["git", "-C", str(config.ROOT_DIR), "rev-parse", "--short", "HEAD"],
            capture_output=True, text=True, timeout=10,
        )
        return out.stdout.strip() or None
    except (OSError, subprocess.SubprocessError):
        return None


def tool_versions(include_python: bool = True) -> Dict[str, Optional[str]]:
    """Versions of the external binaries and key libraries, for the manifest.

    Only reports what is actually importable: the web venv has no torch, and calling this
    from there must not explode.
    """
    versions: Dict[str, Optional[str]] = {
        "hmmer": _cmd_version([config.HMMBUILD_BIN, "-h"], match="HMMER"),
        "mmseqs2": _cmd_version([config.MMSEQS_BIN, "version"]),
    }
    if include_python:
        import platform
        versions["python"] = platform.python_version()
        for mod in ("torch", "transformers", "sklearn", "numpy", "biotite", "gemmi"):
            try:
                versions[mod] = __import__(mod).__version__
            except Exception:
                versions[mod] = None
    return versions


@dataclass
class Manifest:
    stage: str
    run_id: int
    label: str
    model_version: Optional[str] = None
    inputs: Dict[str, str] = field(default_factory=dict)
    outputs: Dict[str, str] = field(default_factory=dict)
    tools: Dict[str, Optional[str]] = field(default_factory=dict)
    n_input: int = 0
    n_output: int = 0
    n_discarded: int = 0
    _t0: float = field(default_factory=time.monotonic)

    def add_input(self, *paths: str | Path) -> None:
        for p in paths:
            self.inputs[str(p)] = sha256_file(p)

    def add_output(self, *paths: str | Path) -> None:
        for p in paths:
            self.outputs[str(p)] = sha256_file(p)

    def counts(self, n_input: int = 0, n_output: int = 0, n_discarded: int = 0) -> None:
        self.n_input, self.n_output, self.n_discarded = n_input, n_output, n_discarded

    @property
    def wall_time_s(self) -> float:
        return time.monotonic() - self._t0

    def write(self) -> Path:
        """Persist to both the manifests table and a JSON file. Returns the file path."""
        payload = {
            "stage": self.stage,
            "run_id": self.run_id,
            "label": self.label,
            "model_version": self.model_version,
            "schema_version": SCHEMA_VERSION,
            "git_commit": git_commit(),
            "input_hashes": self.inputs,
            "output_hashes": self.outputs,
            "tool_versions": self.tools,
            "counts": {
                "n_input": self.n_input,
                "n_output": self.n_output,
                "n_discarded": self.n_discarded,
            },
            "wall_time_s": round(self.wall_time_s, 2),
            "written_at": now(),
        }

        config.MANIFEST_DIR.mkdir(parents=True, exist_ok=True)
        stamp = time.strftime("%Y%m%dT%H%M%SZ", time.gmtime())
        path = config.MANIFEST_DIR / f"{self.stage}_{stamp}.json"
        path.write_text(json.dumps(payload, indent=2, sort_keys=True))

        def _do() -> None:
            with connect() as conn:
                conn.execute(
                    "INSERT INTO manifests(run_id, stage, input_hashes_json, "
                    "output_hashes_json, tool_versions_json, model_version, schema_version, "
                    "git_commit, wall_time_s, written_at) VALUES (?,?,?,?,?,?,?,?,?,?)",
                    (self.run_id, self.stage,
                     json.dumps(self.inputs), json.dumps(self.outputs),
                     json.dumps(self.tools), self.model_version, SCHEMA_VERSION,
                     payload["git_commit"], payload["wall_time_s"], payload["written_at"]),
                )
        retry_write(_do)
        return path


@contextmanager
def stage_manifest(
    stage: str,
    label: str,
    inputs: Optional[Iterable[str | Path]] = None,
    model_version: Optional[str] = None,
    params: Optional[Dict[str, Any]] = None,
) -> Iterator[Manifest]:
    """Open a run, yield a Manifest, and close both on the way out.

    On an exception the run is marked 'error' and the manifest is still written, so a
    failed stage leaves a provenance trail rather than nothing.
    """
    run_id = start_run(stage, label, json.dumps(params) if params else None)
    m = Manifest(stage=stage, run_id=run_id, label=label, model_version=model_version)
    m.tools = tool_versions()
    if inputs:
        m.add_input(*inputs)
    try:
        yield m
    except Exception:
        m.write()
        fail_run(run_id)
        raise
    else:
        m.write()
        finish_run(run_id, m.n_input, m.n_output, m.n_discarded)
