"""Build the profile HMMs that drive the recall stage.

Spec section 5.1: profile HMMs from characterised PET hydrolases plus the relevant ESTHER
families, MMseqs2 for speed across large metagenomic assemblies, HMMER for the final
sensitive pass, recording E-value and profile identity for every candidate. Those
retrieval numbers are the baseline the learned model must beat, so this stage produces
the yardstick as much as it produces candidates.

One profile per positive cluster rather than one profile over everything. The positives
span 11 clusters at 30% identity; a single profile pooling them would be broad and
washed-out, and the per-cluster E-values would no longer mean anything comparable.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .. import config, seqtools

MAFFT_BIN = "mafft"

# A cluster below this size cannot support a meaningful profile: hmmbuild will happily
# build one from two sequences, but its emission probabilities are near-noise.
MIN_CLUSTER_SIZE = 3


class ProfileError(RuntimeError):
    pass


@dataclass
class Profile:
    name: str
    hmm_path: Path
    msa_path: Path
    n_sequences: int
    length: int          # match states in the model
    members: List[str]


def _run(args: List[str], stdout_path: Optional[Path] = None) -> subprocess.CompletedProcess:
    if shutil.which(args[0]) is None:
        raise ProfileError(f"{args[0]} not on PATH")
    if stdout_path is not None:
        with stdout_path.open("w") as fh:
            proc = subprocess.run(args, stdout=fh, stderr=subprocess.PIPE, text=True)
    else:
        proc = subprocess.run(args, capture_output=True, text=True)
    if proc.returncode != 0:
        tail = (proc.stderr or "").strip().splitlines()[-12:]
        raise ProfileError(f"{args[0]} failed ({proc.returncode}):\n" + "\n".join(tail))
    return proc


def align(records: Sequence[Tuple[str, str]], out_msa: Path) -> Path:
    """MAFFT alignment. --auto picks the algorithm from the input size."""
    fasta = out_msa.with_suffix(".unaligned.fasta")
    seqtools.write_fasta(records, fasta)
    out_msa.parent.mkdir(parents=True, exist_ok=True)
    _run([MAFFT_BIN, "--auto", "--quiet", "--anysymbol", str(fasta)], stdout_path=out_msa)
    return out_msa


def hmmbuild(msa: Path, out_hmm: Path, name: str) -> Path:
    out_hmm.parent.mkdir(parents=True, exist_ok=True)
    _run([config.HMMBUILD_BIN, "--amino", "-n", name, str(out_hmm), str(msa)])
    return out_hmm


def _model_length(hmm_path: Path) -> int:
    for line in hmm_path.read_text().splitlines():
        if line.startswith("LENG"):
            return int(line.split()[1])
    return 0


def build_from_clusters(records: Dict[str, str], clusters: Dict[str, str],
                        out_dir: Path, prefix: str = "PLC",
                        min_size: int = MIN_CLUSTER_SIZE) -> Tuple[List[Profile], Dict[str, int]]:
    """Build one profile per cluster.

    `records` maps sequence id to sequence; `clusters` maps member id to representative id
    (the shape mmseqs easy-cluster returns).
    """
    grouped: Dict[str, List[str]] = {}
    for member, rep in clusters.items():
        grouped.setdefault(rep, []).append(member)

    out_dir.mkdir(parents=True, exist_ok=True)
    profiles: List[Profile] = []
    skipped = {"too_small": 0, "n_skipped_members": 0}

    for i, (rep, members) in enumerate(
            sorted(grouped.items(), key=lambda kv: -len(kv[1])), start=1):
        members = [m for m in members if m in records]
        if len(members) < min_size:
            skipped["too_small"] += 1
            skipped["n_skipped_members"] += len(members)
            continue
        name = f"{prefix}_{i:02d}"
        msa = align([(m, records[m]) for m in members], out_dir / f"{name}.afa")
        hmm = hmmbuild(msa, out_dir / f"{name}.hmm", name)
        profiles.append(Profile(name=name, hmm_path=hmm, msa_path=msa,
                                n_sequences=len(members), length=_model_length(hmm),
                                members=members))
    return profiles, skipped


def press(profiles: Sequence[Profile], out_db: Path) -> Path:
    """Concatenate the profiles into one searchable library and hmmpress it.

    hmmpress refuses to overwrite existing binaries, so stale ones are cleared first:
    otherwise a rebuilt library silently keeps searching against the previous models.
    """
    out_db.parent.mkdir(parents=True, exist_ok=True)
    out_db.write_text("".join(p.hmm_path.read_text() for p in profiles))
    for suffix in (".h3f", ".h3i", ".h3m", ".h3p"):
        Path(str(out_db) + suffix).unlink(missing_ok=True)
    _run(["hmmpress", str(out_db)])
    return out_db
