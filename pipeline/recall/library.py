"""The profile library: one HMM per sequence cluster, each with its own catalytic anchor.

This replaces the single pooled profile used for the Phase 2a validation. That profile
proved the triad-reading mechanism works, but it scored 0/111 on the near misses, not
because classic cutinases lack a triad (they plainly have one) but because they never
aligned to a Polyesterase-lipase-cutinase profile well enough for the columns to map.
Since the near misses exist to define the decision boundary (spec section 5.2), a filter
that silently discards them removes exactly what the model most needs to see.

The design here:

  1. Cluster positives AND near misses together at 30% identity. Near misses get their own
     clusters and therefore their own profiles, instead of being forced through a profile
     built from a family they do not belong to.
  2. Give every cluster an anchor: a member whose catalytic triad positions come from
     UniProt's Active site annotation (see anchors.py). A cluster with no annotated member
     gets no profile, and that is reported rather than silently skipped.
  3. Score each candidate against the whole library with hmmscan, take its best-scoring
     profile, and read the triad from THAT profile's anchor.

So a candidate is judged against the family it actually resembles, and the reported
E-value comes from the same profile that supplied its triad call.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .. import config, seqtools
from . import anchors, profiles, triad

HMMSCAN_BIN = "hmmscan"


@dataclass
class LibraryEntry:
    name: str
    profile: profiles.Profile
    anchor: anchors.Anchor


@dataclass
class Library:
    entries: Dict[str, LibraryEntry] = field(default_factory=dict)
    db_path: Optional[Path] = None
    skipped: Dict[str, List[str]] = field(default_factory=dict)

    def __len__(self) -> int:
        return len(self.entries)


def build(sequences: Dict[str, str], accessions: Dict[str, Optional[str]],
          out_dir: Path, identity: float = 0.3, min_cluster_size: int = 3,
          priority: Optional[Sequence[str]] = None) -> Library:
    """Cluster, build a profile per cluster, and anchor each one.

    `priority` lists sequence ids to prefer when choosing a cluster's anchor (the curated
    wild types, whose annotations are the most trustworthy).
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    priority = list(priority or [])

    fasta = seqtools.write_fasta(sequences.items(), out_dir / "all.fasta")
    clusters = seqtools.cluster(fasta, min_seq_id=identity)

    grouped: Dict[str, List[str]] = {}
    for member, rep in clusters.items():
        grouped.setdefault(rep, []).append(member)

    lib = Library(skipped={"too_small": [], "no_anchor": []})

    for i, (rep, members) in enumerate(
            sorted(grouped.items(), key=lambda kv: -len(kv[1])), start=1):
        members = [m for m in members if m in sequences]
        if len(members) < min_cluster_size:
            lib.skipped["too_small"].append(f"{rep} ({len(members)})")
            continue

        # Curated entries first: their UniProt annotations are the ones we trust most.
        ordered = ([m for m in members if m in priority]
                   + [m for m in members if m not in priority])
        anchor = anchors.find_anchor_for_cluster(ordered, sequences, accessions)
        if anchor is None:
            lib.skipped["no_anchor"].append(f"{rep} ({len(members)} members)")
            continue

        name = f"PROF_{i:02d}"
        msa = profiles.align([(m, sequences[m]) for m in members], out_dir / f"{name}.afa")
        hmm = profiles.hmmbuild(msa, out_dir / f"{name}.hmm", name)
        prof = profiles.Profile(name=name, hmm_path=hmm, msa_path=msa,
                                n_sequences=len(members),
                                length=profiles._model_length(hmm), members=members)
        lib.entries[name] = LibraryEntry(name=name, profile=prof, anchor=anchor)

    if lib.entries:
        lib.db_path = profiles.press([e.profile for e in lib.entries.values()],
                                     out_dir / "library.hmm")
    return lib


def hmmscan_best(lib: Library, records: Sequence[Tuple[str, str]],
                 work_dir: Optional[Path] = None,
                 evalue_cutoff: float = 10.0) -> Dict[str, Tuple[str, float, float]]:
    """Best-scoring profile per sequence: {seq_id: (profile_name, evalue, bitscore)}.

    Uses the --tblout table rather than the human-readable output, which is formatted for
    reading and truncates long names.
    """
    if lib.db_path is None:
        return {}
    # An assembly where the prefilter found nothing is a normal outcome, not an error:
    # several gut assemblies contain no polyesterase-like protein at all. hmmscan treats
    # an empty input file as malformed and exits non-zero, which crashed a 50-file scan
    # on file 46. Nothing to scan means no hits, so say so and return.
    if not records:
        return {}
    if shutil.which(HMMSCAN_BIN) is None:
        raise profiles.ProfileError(f"{HMMSCAN_BIN} not on PATH")

    work_dir = work_dir or config.INTERIM_DIR / "scan"
    work_dir.mkdir(parents=True, exist_ok=True)
    fasta = seqtools.write_fasta(records, work_dir / "scan_in.fasta")
    tbl = work_dir / "scan.tbl"

    proc = subprocess.run(
        [HMMSCAN_BIN, "--tblout", str(tbl), "-E", str(evalue_cutoff),
         "--noali", str(lib.db_path), str(fasta)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise profiles.ProfileError(
            f"hmmscan failed: {(proc.stderr or '').strip()[:400]}")

    best: Dict[str, Tuple[str, float, float]] = {}
    for line in tbl.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split()
        if len(parts) < 6:
            continue
        prof_name, seq_id, evalue, score = parts[0], parts[2], float(parts[4]), float(parts[5])
        if seq_id not in best or score > best[seq_id][2]:
            best[seq_id] = (prof_name, evalue, score)
    return best


def call_triads(lib: Library, records: Sequence[Tuple[str, str]],
                work_dir: Optional[Path] = None
                ) -> Tuple[Dict[str, triad.TriadCall], Dict[str, Tuple[str, float, float]]]:
    """Assign each sequence to its best profile, then read the triad from that anchor.

    Sequences matching no profile are absent from the returned calls: they were not
    judged, which is different from being judged and failing, and the caller must report
    the two separately.
    """
    if not records:
        return {}, {}
    best = hmmscan_best(lib, records, work_dir=work_dir)
    seqs = dict(records)

    by_profile: Dict[str, List[str]] = {}
    for sid, (prof_name, _e, _s) in best.items():
        by_profile.setdefault(prof_name, []).append(sid)

    calls: Dict[str, triad.TriadCall] = {}
    work_dir = work_dir or config.INTERIM_DIR / "scan"
    for prof_name, members in by_profile.items():
        entry = lib.entries.get(prof_name)
        if entry is None:
            continue
        anchor = entry.anchor
        # The anchor must be in the alignment for its columns to be readable, so it is
        # added explicitly rather than assumed to be among the candidates.
        recs = [(anchor.sequence_id, anchor.sequence)] + [
            (m, seqs[m]) for m in members if m != anchor.sequence_id
        ]
        sub = triad.call_triads(
            entry.profile.hmm_path, recs, reference_id=anchor.sequence_id,
            work_dir=work_dir / prof_name,
            triad={"ser": anchor.ser, "asp": anchor.asp, "his": anchor.his},
        )
        for m in members:
            if m in sub:
                calls[m] = sub[m]
    return calls, best
