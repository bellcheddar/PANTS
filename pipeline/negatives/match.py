"""Select a matched hard-negative set and write it to `characterised_enzymes`.

Spec section 5.2 requires negatives matched to positives on three axes:

  1. Sequence length distribution: negatives are drawn to follow the positives' own
     length histogram, so length carries no signal.
  2. Identity to the nearest positive: selection PREFERS high-identity negatives, because
     a negative that looks nothing like a positive is not a hard negative and teaches the
     model nothing about the boundary.
  3. Taxonomic breadth: a per-genus cap, so the set cannot be won by recognising one
     over-represented organism.

Everything selected here is a *characterised* member of a family with no polyester
activity, so these are genuine negatives rather than the unlabelled pool. The unlabelled
metagenomic candidates are handled as PU in the training stage, which is a different
problem (spec section 5.3).
"""

from __future__ import annotations

import json
import random
from collections import Counter, defaultdict
from typing import Dict, List, Optional, Tuple

from .. import config, seqtools
from ..db import connect, now
from . import esther

# Default target size. With order 10^1 sequence-resolved positives this is a deliberately
# generous ratio: hard negatives are cheap, and the PU-corrected loss handles imbalance.
DEFAULT_N_NEGATIVES = 400
DEFAULT_GENUS_CAP = 8


def _genus(organism: Optional[str]) -> str:
    return (organism or "unknown").split()[0]


def positives_from_db() -> List[Tuple[str, str]]:
    """Sequence-resolved positives in the PETase/cutinase families.

    MHETase is excluded: it is a Tannase-family enzyme with a lid domain, a different
    problem with its own seed and negatives in v2. Including it here would blur the very
    boundary this negative set exists to define.
    """
    with connect() as conn:
        return [
            (r["enzyme_id"], r["sequence"])
            for r in conn.execute(
                "SELECT enzyme_id, sequence FROM characterised_enzymes "
                "WHERE is_positive=1 AND sequence IS NOT NULL AND family != 'mhetase_like'"
            )
        ]


def select(hits: List[esther.EstherHit], positives: List[Tuple[str, str]],
           n_target: int = DEFAULT_N_NEGATIVES, genus_cap: int = DEFAULT_GENUS_CAP,
           seed: int = 0) -> Tuple[List[esther.EstherHit], Dict[str, object]]:
    """Choose a length-matched, identity-preferring, taxonomically capped subset."""
    rng = random.Random(seed)

    # --- identity of every candidate negative to its nearest positive ---
    tmp = config.INTERIM_DIR
    tmp.mkdir(parents=True, exist_ok=True)
    neg_fa = seqtools.write_fasta(((h.accession, h.sequence) for h in hits), tmp / "neg_all.fasta")
    pos_fa = seqtools.write_fasta(positives, tmp / "pos.fasta")
    best = seqtools.best_identity_to(neg_fa, pos_fa)   # {acc: (positive_id, fident)}

    identity = {h.accession: best.get(h.accession, ("", 0.0))[1] for h in hits}
    nearest = {h.accession: best.get(h.accession, ("", 0.0))[0] for h in hits}

    # --- length matching: follow the positives' own histogram ---
    pos_lengths = [len(s) for _, s in positives]
    by_length: Dict[int, List[esther.EstherHit]] = defaultdict(list)
    for h in hits:
        by_length[h.length].append(h)
    # Prefer high identity within each length bucket: those are the hard ones.
    for bucket in by_length.values():
        bucket.sort(key=lambda h: identity[h.accession], reverse=True)

    chosen: List[esther.EstherHit] = []
    taken = set()
    genus_count: Counter = Counter()
    per_positive = max(1, n_target // max(1, len(pos_lengths)))

    for target_len in pos_lengths:
        picked = 0
        # Walk outwards from the positive's exact length, so the resulting distribution is
        # centred on it rather than on the negatives' own (much broader) distribution.
        for delta in range(0, 60):
            for cand_len in {target_len - delta, target_len + delta}:
                for h in by_length.get(cand_len, []):
                    if picked >= per_positive:
                        break
                    if h.accession in taken:
                        continue
                    g = _genus(h.organism)
                    if genus_count[g] >= genus_cap:
                        continue
                    chosen.append(h)
                    taken.add(h.accession)
                    genus_count[g] += 1
                    picked += 1
            if picked >= per_positive:
                break

    rng.shuffle(chosen)

    chosen_lengths = [h.length for h in chosen]
    chosen_ident = [identity[h.accession] for h in chosen]
    report = {
        "n_chosen": len(chosen),
        "n_pool": len(hits),
        "positive_length_range": (min(pos_lengths), max(pos_lengths)),
        "chosen_length_range": (min(chosen_lengths), max(chosen_lengths)) if chosen else None,
        "chosen_length_mean": round(sum(chosen_lengths) / len(chosen), 1) if chosen else None,
        "positive_length_mean": round(sum(pos_lengths) / len(pos_lengths), 1),
        "identity_mean": round(sum(chosen_ident) / len(chosen), 4) if chosen else None,
        "identity_max": round(max(chosen_ident), 4) if chosen else None,
        "n_with_any_identity": sum(1 for v in chosen_ident if v > 0),
        "n_genera": len({_genus(h.organism) for h in chosen}),
        "top_genera": Counter(_genus(h.organism) for h in chosen).most_common(8),
        "family_breakdown": Counter(h.family for h in chosen).most_common(),
    }
    return chosen, {"identity": identity, "nearest": nearest, "report": report}


def write(chosen: List[esther.EstherHit], near_misses: List[esther.EstherHit],
          identity: Dict[str, float], nearest: Dict[str, str]) -> int:
    """Upsert negatives and near misses into characterised_enzymes."""
    n = 0
    with connect() as conn:
        for role, group in (("negative", chosen), ("near_miss", near_misses)):
            for h in group:
                conn.execute(
                    "INSERT INTO characterised_enzymes "
                    "(enzyme_id, uniprot, organism, family, sequence, seq_length, "
                    " is_positive, is_negative, is_near_miss, esther_family, "
                    " matched_positive_id, activity_substrate_notes, source_ref, added_at) "
                    "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?) "
                    "ON CONFLICT(enzyme_id) DO UPDATE SET "
                    " sequence=excluded.sequence, seq_length=excluded.seq_length, "
                    " is_negative=excluded.is_negative, is_near_miss=excluded.is_near_miss, "
                    " esther_family=excluded.esther_family, "
                    " matched_positive_id=excluded.matched_positive_id",
                    (
                        f"ESTHER:{h.accession}", h.accession, h.organism,
                        "cutinase" if role == "near_miss" else "other",
                        h.sequence, h.length,
                        0, int(role == "negative"), int(role == "near_miss"),
                        h.family, nearest.get(h.accession) or None,
                        (f"ESTHER family {h.family}. "
                         f"Identity to nearest positive: {identity.get(h.accession, 0.0):.3f}. "
                         + ("Near miss: active on soluble esters, not a meaningful "
                            "degrader of crystalline PET." if role == "near_miss"
                            else "No polyester activity.")),
                        "ESTHER/UniProt", now(),
                    ),
                )
                n += 1
    return n
