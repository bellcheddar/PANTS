"""The recall stage end to end: metagenome FASTA in, triad-complete candidates out.

Spec section 5.1's two-stage design, and the reason for it: MMseqs2 casts the net across
millions of sequences fast, HMMER makes the sensitive final call on the survivors. Running
HMMER over the raw metagenome would work and would take far longer for the same answer.

Every candidate keeps its retrieval numbers (E-value, bitscore, profile identity), because
those are the baseline the learned model has to beat (spec section 8). The recall stage
therefore produces the yardstick as well as the candidates.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

from .. import config, seqtools
from ..db import connect, now
from ..db.manifest import stage_manifest
from . import library as lib_mod

STAGE = "recall"


def candidate_id(sequence: str) -> str:
    """Stable id: 'PANTS-' plus the first 12 hex of SHA256(sequence).

    Content-addressed on purpose, so the same protein recovered from a different assembly
    collapses to one row instead of becoming a duplicate candidate.
    """
    return "PANTS-" + hashlib.sha256(sequence.encode()).hexdigest()[:12]


def iter_fasta(path: Path) -> Iterable[Tuple[str, str]]:
    """Stream (id, sequence) without loading a multi-GB file into memory."""
    rid: Optional[str] = None
    buf: List[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.rstrip()
            if line.startswith(">"):
                if rid is not None:
                    yield rid, "".join(buf)
                rid, buf = line[1:].split()[0], []
            elif rid is not None:
                buf.append(line)
    if rid is not None:
        yield rid, "".join(buf)


def prefilter(metagenome_fasta: Path, positives: List[Tuple[str, str]],
              work_dir: Path, sensitivity: Optional[float] = None,
              max_evalue: float = 1e-5) -> Dict[str, Tuple[str, float, float]]:
    """MMseqs2 pass: metagenome as query, positives as target.

    This direction matters. Querying the metagenome against a small target database gives
    one best-hit row per metagenomic protein, which is what the candidate table needs;
    the reverse direction returns the positives' best hits and silently caps how many
    metagenomic proteins can ever be recovered at --max-seqs.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    pos_fa = seqtools.write_fasta(positives, work_dir / "positives.fasta")
    hits = seqtools.easy_search(metagenome_fasta, pos_fa, sensitivity=sensitivity)

    best: Dict[str, Tuple[str, float, float]] = {}
    for h in hits:
        q, t, ev, fid = str(h["query"]), str(h["target"]), float(h["evalue"]), float(h["fident"])
        if ev > max_evalue:
            continue
        if q not in best or ev < best[q][1]:
            best[q] = (t, ev, fid)
    return best


@dataclass
class RecallResult:
    n_scanned: int
    n_prefilter: int
    n_profile_matched: int
    n_triad_complete: int
    n_written: int
    discarded: Dict[str, int]


def run(metagenome_fastas: List[Path], lib: lib_mod.Library,
        positives: List[Tuple[str, str]], environment: str,
        label: str = "v1", max_evalue: float = 1e-5,
        work_dir: Optional[Path] = None) -> RecallResult:
    """Full recall over one or more metagenome FASTA files."""
    work_dir = work_dir or config.INTERIM_DIR / "recall"
    work_dir.mkdir(parents=True, exist_ok=True)
    discarded = {"prefilter": 0, "no_profile": 0, "incomplete_triad": 0}

    with stage_manifest(STAGE, label=label, inputs=metagenome_fastas,
                        params={"max_evalue": max_evalue, "environment": environment}) as m:
        n_scanned = 0
        surviving: Dict[str, str] = {}
        retrieval: Dict[str, Tuple[str, float, float]] = {}

        for fasta in metagenome_fastas:
            seqs = dict(iter_fasta(fasta))
            n_scanned += len(seqs)
            best = prefilter(fasta, positives, work_dir, max_evalue=max_evalue)
            discarded["prefilter"] += len(seqs) - len(best)
            for sid, hit in best.items():
                if sid in seqs:
                    surviving[f"{fasta.stem}|{sid}"] = seqs[sid]
                    retrieval[f"{fasta.stem}|{sid}"] = hit

        n_prefilter = len(surviving)

        # Sensitive pass: assign each survivor to its best profile and read the triad
        # from that profile's own anchor.
        calls, scan = lib_mod.call_triads(lib, list(surviving.items()),
                                          work_dir=work_dir / "scan")
        discarded["no_profile"] = n_prefilter - len(calls)
        complete = {k: v for k, v in calls.items() if v.complete}
        discarded["incomplete_triad"] = len(calls) - len(complete)

        n_written = 0
        with connect() as conn:
            for key, call in complete.items():
                seq = surviving[key]
                cid = candidate_id(seq)
                target, evalue, fident = retrieval.get(key, ("", None, None))
                prof, prof_e, prof_score = scan.get(key, (None, None, None))
                contig = key.split("|", 1)[1]
                conn.execute(
                    "INSERT INTO candidates (candidate_id, source_environment, assembly_id, "
                    " contig_id, sequence, seq_length, has_complete_triad, "
                    " triad_positions_json, recall_method, recall_evalue, recall_bitscore, "
                    " recall_profile_identity, nearest_characterised_id, "
                    " nearest_characterised_identity, first_seen_at) "
                    "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) "
                    "ON CONFLICT(candidate_id) DO NOTHING",
                    (cid, environment, key.split("|", 1)[0], contig, seq, len(seq), 1,
                     f'{{"ser": {call.ser_pos}, "asp": {call.asp_pos}, "his": {call.his_pos}}}',
                     "mmseqs2+hmmscan", prof_e, prof_score, fident, target, fident, now()),
                )
                n_written += 1

        m.counts(n_input=n_scanned, n_output=n_written,
                 n_discarded=n_scanned - n_written)

    return RecallResult(n_scanned=n_scanned, n_prefilter=n_prefilter,
                        n_profile_matched=len(calls), n_triad_complete=len(complete),
                        n_written=n_written, discarded=discarded)
