"""Catalytic triad and oxyanion hole detection, by profile alignment.

Spec section 5.1 gates the whole pipeline on this: "Filter for a complete catalytic triad
and a recognisable oxyanion hole before scoring. Report how many candidates are discarded
at this step."

The method, and why it is not a motif regex: the triad is Ser-Asp-His in sequence order
but the three residues sit far apart in sequence and are only brought together by the
fold. A regex over GxSxG finds the nucleophile elbow and nothing else, and would pass any
α/β-hydrolase while saying nothing about whether the other two partners are present and
correctly placed.

Instead every candidate is aligned to the profile HMM with `hmmalign`, and the columns
corresponding to IsPETase's OWN verified triad positions are read off. IsPETase's
S160/D206/H237 were checked directly against its UniProt sequence during curation, so the
reference is anchored to something confirmed rather than assumed. A candidate has a
complete triad when it presents Ser, Asp and His at those three profile columns.

That also handles insertions and deletions correctly: alignment columns are the shared
coordinate system, so a candidate with a 12-residue insertion before the catalytic Ser
still maps onto the right column, where residue-offset arithmetic would silently drift.
"""

from __future__ import annotations

import shutil
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

from .. import config

HMMALIGN_BIN = "hmmalign"

# The catalytic triad in IsPETase (UniProt A0A0K8P6T7) precursor numbering. Verified
# against the fetched sequence during Phase 1 curation: position 160 reads S, 206 reads D,
# 237 reads H. Do not change these without re-running that check.
ISPETASE_TRIAD = {"ser": 160, "asp": 206, "his": 237}

# The oxyanion hole in the PETase fold is formed by backbone amides of the residue
# following the nucleophile elbow and of a residue near the conserved Trp. Backbone
# geometry is a structural property, so sequence alone can only support a WEAK proxy:
# we record the residues at those columns rather than claiming to detect the hole.
# The real determination happens in the structure stage (spec section 6.3).
ISPETASE_OXYANION = {"y87": 87, "m161": 161}


class AlignError(RuntimeError):
    pass


@dataclass
class TriadCall:
    sequence_id: str
    ser: Optional[str]
    asp: Optional[str]
    his: Optional[str]
    ser_pos: Optional[int]
    asp_pos: Optional[int]
    his_pos: Optional[int]
    oxyanion_residues: Dict[str, Optional[str]]
    aligned_fraction: float

    @property
    def complete(self) -> bool:
        return (self.ser, self.asp, self.his) == ("S", "D", "H")

    @property
    def reason(self) -> str:
        if self.complete:
            return "complete"
        missing = [n for n, v, want in
                   (("Ser", self.ser, "S"), ("Asp", self.asp, "D"), ("His", self.his, "H"))
                   if v != want]
        return "missing/substituted: " + ", ".join(
            f"{n}={'gap' if dict(Ser=self.ser, Asp=self.asp, His=self.his)[n] is None else dict(Ser=self.ser, Asp=self.asp, His=self.his)[n]}"
            for n in missing)


def _read_stockholm(path: Path) -> Dict[str, str]:
    """Parse hmmalign's Stockholm output into {id: aligned_sequence}.

    Stockholm interleaves blocks for long alignments, so sequence lines for one id must be
    concatenated across blocks rather than overwritten.
    """
    seqs: Dict[str, List[str]] = {}
    for line in path.read_text().splitlines():
        line = line.rstrip()
        if not line or line.startswith("#") or line == "//":
            continue
        parts = line.split(None, 1)
        if len(parts) != 2:
            continue
        seqs.setdefault(parts[0], []).append(parts[1].strip())
    return {k: "".join(v) for k, v in seqs.items()}


def hmmalign(hmm: Path, fasta: Path, out_sto: Path) -> Dict[str, str]:
    if shutil.which(HMMALIGN_BIN) is None:
        raise AlignError(f"{HMMALIGN_BIN} not on PATH")
    out_sto.parent.mkdir(parents=True, exist_ok=True)
    proc = subprocess.run(
        [HMMALIGN_BIN, "--amino", "-o", str(out_sto), str(hmm), str(fasta)],
        capture_output=True, text=True,
    )
    if proc.returncode != 0:
        raise AlignError(f"hmmalign failed: {(proc.stderr or '').strip()[:400]}")
    return _read_stockholm(out_sto)


def _match_columns(aligned: str) -> List[int]:
    """Indices of alignment columns that are MATCH states.

    hmmalign writes match/delete states in upper case and dots, and insertions relative to
    the model in lower case. Only the upper-case/dash columns form the model's coordinate
    system shared across all sequences.
    """
    return [i for i, c in enumerate(aligned) if c.isupper() or c == "-"]


def _residue_at_model_position(aligned: str, ref_aligned: str, ref_pos: int
                               ) -> Tuple[Optional[str], Optional[int]]:
    """Residue in `aligned` at the column where `ref_aligned` holds its `ref_pos`-th residue.

    Returns (residue or None if gapped, 1-based position in the ungapped candidate).
    """
    # Count EVERY alphabetic character, upper and lower. hmmalign writes residues that
    # are insertions relative to the model in lower case, but they are still residues of
    # the sequence and still consume a position in its own numbering. Counting only the
    # upper-case (match-state) columns walks off by the number of insertions before the
    # target, which is exactly what made IsPETase report its catalytic Ser160 as A171.
    seen = 0
    col = None
    for i, c in enumerate(ref_aligned):
        if c.isalpha():
            seen += 1
            if seen == ref_pos:
                col = i
                break
    if col is None:
        return None, None

    ch = aligned[col] if col < len(aligned) else None
    if ch is None or not ch.isalpha():
        return None, None
    pos = sum(1 for c in aligned[:col + 1] if c.isalpha())
    return ch.upper(), pos


def call_triads(hmm: Path, records: Sequence[Tuple[str, str]], reference_id: str,
                work_dir: Optional[Path] = None,
                triad: Dict[str, int] = None,
                oxyanion: Dict[str, int] = None) -> Dict[str, TriadCall]:
    """Align every record to the profile and read the triad off the reference's columns.

    `reference_id` must be present in `records` and must be the sequence whose numbering
    `triad` refers to (IsPETase by default).
    """
    from .. import seqtools

    triad = triad or ISPETASE_TRIAD
    oxyanion = oxyanion or ISPETASE_OXYANION
    work_dir = work_dir or config.INTERIM_DIR / "triad"
    work_dir.mkdir(parents=True, exist_ok=True)

    ids = [r[0] for r in records]
    if reference_id not in ids:
        raise AlignError(f"reference {reference_id!r} not among the sequences to align")

    fasta = seqtools.write_fasta(records, work_dir / "to_align.fasta")
    aligned = hmmalign(hmm, fasta, work_dir / "aligned.sto")

    ref = aligned.get(reference_id)
    if ref is None:
        raise AlignError(f"reference {reference_id!r} missing from the alignment")

    calls: Dict[str, TriadCall] = {}
    for sid, _ in records:
        al = aligned.get(sid)
        if al is None:
            continue
        ser, ser_p = _residue_at_model_position(al, ref, triad["ser"])
        asp, asp_p = _residue_at_model_position(al, ref, triad["asp"])
        his, his_p = _residue_at_model_position(al, ref, triad["his"])
        oxy = {name: _residue_at_model_position(al, ref, pos)[0]
               for name, pos in oxyanion.items()}
        n_res = sum(1 for c in al if c.isalpha())
        n_match = sum(1 for c in al if c.isupper())
        calls[sid] = TriadCall(
            sequence_id=sid, ser=ser, asp=asp, his=his,
            ser_pos=ser_p, asp_pos=asp_p, his_pos=his_p,
            oxyanion_residues=oxy,
            aligned_fraction=round(n_match / n_res, 3) if n_res else 0.0,
        )
    return calls
