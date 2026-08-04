"""Thin wrappers over the MMseqs2 and HMMER binaries.

Shelled out to rather than bound as libraries, matching BoltzMaker's house style for
external tools: the CLI is the stable interface, and the version that ran is recorded in
the manifest.
"""

from __future__ import annotations

import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from . import config


class ToolError(RuntimeError):
    """An external tool exited non-zero."""


def _run(args: List[str]) -> subprocess.CompletedProcess:
    proc = subprocess.run(args, capture_output=True, text=True)
    if proc.returncode != 0:
        tail = (proc.stderr or proc.stdout or "").strip().splitlines()[-15:]
        raise ToolError(f"{args[0]} failed ({proc.returncode}):\n" + "\n".join(tail))
    return proc


def write_fasta(records: Iterable[Tuple[str, str]], path: str | Path, width: int = 60) -> Path:
    """Write (id, sequence) pairs to FASTA."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n")
            for i in range(0, len(seq), width):
                fh.write(seq[i:i + width] + "\n")
    return path


def read_fasta(path: str | Path) -> Dict[str, str]:
    seqs: Dict[str, str] = {}
    rid: Optional[str] = None
    buf: List[str] = []
    with Path(path).open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if rid is not None:
                    seqs[rid] = "".join(buf)
                rid, buf = line[1:].split()[0], []
            elif rid is not None:
                buf.append(line)
    if rid is not None:
        seqs[rid] = "".join(buf)
    return seqs


def easy_search(query_fasta: str | Path, target_fasta: str | Path,
                sensitivity: Optional[float] = None,
                max_seqs: int = 1000) -> List[Dict[str, object]]:
    """MMseqs2 easy-search. Returns hit dicts with query, target, fident, evalue, bits.

    fident is a FRACTION (0 to 1), not a percentage: MMseqs2 reports it that way and
    treating it as a percentage silently makes every identity filter a no-op.
    """
    sensitivity = config.MMSEQS_SENSITIVITY if sensitivity is None else sensitivity
    if shutil.which(config.MMSEQS_BIN) is None:
        raise ToolError(f"{config.MMSEQS_BIN} not on PATH")

    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "hits.m8"
        _run([
            config.MMSEQS_BIN, "easy-search", str(query_fasta), str(target_fasta),
            str(out), str(Path(tmp) / "mmseqs_tmp"),
            "-s", str(sensitivity), "--max-seqs", str(max_seqs),
            "--format-output", "query,target,fident,alnlen,evalue,bits",
            "-v", "1",
        ])
        hits: List[Dict[str, object]] = []
        if out.exists():
            for line in out.read_text().splitlines():
                if not line.strip():
                    continue
                q, t, fident, alnlen, evalue, bits = line.split("\t")
                hits.append({
                    "query": q, "target": t, "fident": float(fident),
                    "alnlen": int(alnlen), "evalue": float(evalue), "bits": float(bits),
                })
        return hits


def best_identity_to(query_fasta: str | Path, target_fasta: str | Path) -> Dict[str, Tuple[str, float]]:
    """For each query, its best target and identity fraction. Queries with no hit are absent."""
    best: Dict[str, Tuple[str, float]] = {}
    for h in easy_search(query_fasta, target_fasta):
        q, t, f = str(h["query"]), str(h["target"]), float(h["fident"])
        if q not in best or f > best[q][1]:
            best[q] = (t, f)
    return best


def cluster(fasta: str | Path, min_seq_id: float, coverage: float = 0.8) -> Dict[str, str]:
    """MMseqs2 easy-cluster at a given identity. Returns {member_id: representative_id}.

    Used for the evaluation splits, which are by CLUSTER and never by sequence
    (spec section 8).
    """
    if shutil.which(config.MMSEQS_BIN) is None:
        raise ToolError(f"{config.MMSEQS_BIN} not on PATH")
    with tempfile.TemporaryDirectory() as tmp:
        prefix = Path(tmp) / "clu"
        _run([
            config.MMSEQS_BIN, "easy-cluster", str(fasta), str(prefix),
            str(Path(tmp) / "mmseqs_tmp"),
            "--min-seq-id", str(min_seq_id), "-c", str(coverage), "-v", "1",
        ])
        mapping: Dict[str, str] = {}
        tsv = Path(f"{prefix}_cluster.tsv")
        if tsv.exists():
            for line in tsv.read_text().splitlines():
                if line.strip():
                    rep, member = line.split("\t")[:2]
                    mapping[member] = rep
        return mapping
