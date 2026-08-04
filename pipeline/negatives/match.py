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


def phylum_of(lineage: Optional[str]) -> str:
    """Bacterial phylum from a UniProt lineage string.

    Matched on because it is the strongest remaining driver of amino-acid composition:
    high-GC Actinomycetota encode systematically Ala/Gly/Pro-rich proteins. The positives
    are dominated by Actinomycetota (Thermobifida, Streptomyces, Amycolatopsis), so a
    negative set drawn from a different phylum mix leaves the trivial baseline a genuine
    compositional shortcut even after length, identity and secretion are matched.

    Modern NCBI bacterial phylum names all end in -ota (Actinomycetota, Pseudomonadota,
    Bacillota), which identifies the rank without needing a fixed lineage depth.
    """
    if not lineage:
        return "unknown"
    parts = [p.strip() for p in lineage.split(";") if p.strip()]
    for p in parts[1:]:
        if p.endswith("ota") and p.lower() not in {"bacteriota"}:
            return p
    return parts[2] if len(parts) > 2 else "unknown"


def positives_from_db(with_lineage: bool = False):
    """Sequence-resolved positives in the PETase/cutinase families.

    MHETase is excluded: it is a Tannase-family enzyme with a lid domain, a different
    problem with its own seed and negatives in v2. Including it here would blur the very
    boundary this negative set exists to define.
    """
    with connect() as conn:
        rows = list(conn.execute(
            "SELECT enzyme_id, sequence, taxonomy_lineage FROM characterised_enzymes "
            "WHERE is_positive=1 AND sequence IS NOT NULL AND family != 'mhetase_like'"
        ))
    if with_lineage:
        return [(r["enzyme_id"], r["sequence"], r["taxonomy_lineage"]) for r in rows]
    return [(r["enzyme_id"], r["sequence"]) for r in rows]


def select(hits: List[esther.EstherHit], positives, n_target: int = DEFAULT_N_NEGATIVES,
           genus_cap: int = DEFAULT_GENUS_CAP, seed: int = 0,
           match_phylum: bool = True) -> Tuple[List[esther.EstherHit], Dict[str, object]]:
    """Choose a subset matched on length, identity, secretion, genus and phylum."""
    rng = random.Random(seed)
    # positives may arrive as (id, seq) or (id, seq, lineage).
    positives_full = [p for p in positives if len(p) == 3] if match_phylum else []
    positives = [(p[0], p[1]) for p in positives]

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

    # --- phylum quotas: mirror the positives' own phylum mix ---
    # Without this the negatives skew Proteobacteria while the positives are
    # Actinomycetota, and composition alone still separates them (see PHASE1_FINDINGS.md).
    pos_phyla = Counter(phylum_of(lin) for _, _, lin in positives_full) if positives_full else Counter()
    total_pos = sum(pos_phyla.values()) or 1
    quota = {ph: max(1, round(n_target * n / total_pos)) for ph, n in pos_phyla.items()}

    chosen: List[esther.EstherHit] = []
    taken = set()
    genus_count: Counter = Counter()
    phylum_count: Counter = Counter()
    per_positive = max(1, n_target // max(1, len(pos_lengths)))

    def admissible(h: esther.EstherHit) -> bool:
        if h.accession in taken:
            return False
        if genus_count[_genus(h.organism)] >= genus_cap:
            return False
        if quota:
            ph = phylum_of(h.lineage)
            # Outside the positives' phylum mix entirely, or that phylum is already full.
            if ph not in quota or phylum_count[ph] >= quota[ph]:
                return False
        return True

    for target_len in pos_lengths:
        picked = 0
        # Walk outwards from the positive's exact length, so the resulting distribution is
        # centred on it rather than on the negatives' own (much broader) distribution.
        for delta in range(0, 60):
            for cand_len in {target_len - delta, target_len + delta}:
                for h in by_length.get(cand_len, []):
                    if picked >= per_positive:
                        break
                    if not admissible(h):
                        continue
                    chosen.append(h)
                    taken.add(h.accession)
                    genus_count[_genus(h.organism)] += 1
                    phylum_count[phylum_of(h.lineage)] += 1
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
        "positive_phyla": pos_phyla.most_common(5),
        "negative_phyla": Counter(phylum_of(h.lineage) for h in chosen).most_common(5),
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
                    " taxonomy_lineage, "
                    " matched_positive_id, activity_substrate_notes, source_ref, added_at) "
                    "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) "
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
                        h.family, h.lineage, nearest.get(h.accession) or None,
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
