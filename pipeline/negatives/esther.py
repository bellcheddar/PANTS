"""ESTHER hard negative harvest.

Spec section 5.2 calls this the single most important dataset decision, and PLAN_v1.md
makes it risk 1: a badly matched negative set lets the model win on a shortcut (length,
taxonomic origin, anything but the discriminative signal) and every downstream number
then looks fine while meaning nothing.

Family assignment comes from the ESTHER cross-reference carried on UniProt entries, read
from the `xref_esther` field's FamilyName property. UniProt has no per-family query
(`xref:esther-<family>` returns nothing), so the harvest streams a length- and
taxonomy-restricted slice and classifies client-side.

Verified empirically during Phase 1, not assumed:
  IsPETase, LCC, TfCut2, Cut190 -> Polyesterase-lipase-cutinase   (the positive family)
  MHETase                       -> Tannase                        (a different family
      entirely, which is exactly why a PETase-seeded profile search never reaches it:
      spec section 2, point 5)
  a SEPARATE `Cutinase` family exists alongside Polyesterase-lipase-cutinase, and is the
      natural home of the near misses the spec asks for by name.
"""

from __future__ import annotations

import re
from collections import Counter
from dataclasses import dataclass
from typing import Dict, Iterator, List, Optional, Set

from .. import config, http

SEARCH_URL = config.UNIPROT_REST_URL
FIELDS = "accession,protein_name,organism_name,organism_id,length,sequence,xref_esther,lineage"
_LINK_NEXT = re.compile(r'<([^>]+)>;\s*rel="next"')

# --------------------------------------------------------------------------------------
# Family policy
# --------------------------------------------------------------------------------------
# The positive family. Never a negative: these are the enzymes we are trying to find.
POSITIVE_FAMILIES: Set[str] = {"Polyesterase-lipase-cutinase"}

# Active on a polyester, just not on PET. Excluded from the negative set entirely rather
# than being called inactive: labelling a PHB depolymerase "no polyester activity" would
# be factually wrong and would teach the model the opposite of what we want.
POLYESTER_ACTIVE_EXCLUDE: Set[str] = {
    "Esterase_phb", "PHB_depolymerase", "PHB_depolymerase_lipase", "Polyhydroxybutyrate_depolymerase",
    "Tannase",   # MHETase's family: a genuine target, handled by its own v2 pipeline
}

# Near misses: fold, triad and esterase chemistry shared, active on soluble esters and on
# cutin, but not meaningful degraders of crystalline PET. Spec section 5.2 asks for these
# explicitly because they define the decision boundary.
NEAR_MISS_FAMILIES: Set[str] = {"Cutinase"}

# Hard negatives: alpha/beta hydrolase fold, Ser-His-Asp triad, no polyester activity.
# An explicit allowlist rather than "everything else", so a family nobody has thought
# about cannot silently become training signal.
NEGATIVE_FAMILIES: Set[str] = {
    # --- bacterial lipases: the classic fold-and-triad sharers ---
    "Bacterial_lip_FamI.1", "Bacterial_lip_FamI.2", "Bacterial_lip_FamI.3",
    "Bacterial_lip_FamI.4", "Bacterial_lip_FamI.5", "Bacterial_lip_FamI.6",
    "Bacterial_lip_FamI.7", "Bacterial_lip_FamI.8",
    # --- carboxylesterases ---
    "Carb_B_Bacteria", "Acetyl_esterase", "BioH", "A85-EsteraseD-FGH",
    "Carboxymethylbutenolide_lactonase",
    # --- hormone-sensitive lipase family (named in spec section 5.2) ---
    "Hormone-sensitive_lipase_like",
    # --- other triad-bearing alpha/beta hydrolases with no polyester activity ---
    "Dienelactone_hydrolase", "Proline_iminopeptidase", "Homoserine_transacetylase",
    "Thioesterase", "Thioesterase_acyl-transferase", "Epoxide_hydrolase",
    "Haloperoxidase", "Carbon-carbon_bond_hydrolase", "MenH_SHCHC", "RutD",
    "AlphaBeta_hydrolase", "6_AlphaBeta_hydrolase", "abh_upf00227",
    "Lactobacillus_peptidase", "Duf_1100-R", "OHBut_olig_hydro_put",
    "A85-Mycolyl-transferase",
}


@dataclass
class EstherHit:
    accession: str
    protein_name: Optional[str]
    organism: Optional[str]
    taxid: Optional[int]
    lineage: Optional[str]
    sequence: str
    length: int
    family: Optional[str]

    @property
    def is_clean(self) -> bool:
        return bool(self.sequence) and set(self.sequence) <= set("ACDEFGHIKLMNPQRSTVWY")

    @property
    def role(self) -> Optional[str]:
        """'negative', 'near_miss', or None if the family is not usable."""
        if self.family in NEGATIVE_FAMILIES:
            return "negative"
        if self.family in NEAR_MISS_FAMILIES:
            return "near_miss"
        return None


def _family_of(result: dict) -> Optional[str]:
    for x in result.get("uniProtKBCrossReferences", []):
        if x.get("database") == "ESTHER":
            for p in x.get("properties", []):
                if p.get("key") == "FamilyName":
                    return p.get("value")
    return None


def stream(query: str, max_results: Optional[int] = None) -> Iterator[EstherHit]:
    """Stream ESTHER-annotated UniProt entries for a query, following Link pagination."""
    url: Optional[str] = SEARCH_URL
    params: Optional[dict] = {"query": query, "format": "json", "size": 500, "fields": FIELDS}
    seen = 0
    while url:
        resp = http.get(url, params=params)
        resp.raise_for_status()
        for res in resp.json().get("results", []):
            org = res.get("organism") or {}
            lineage = org.get("lineage")
            seq = (res.get("sequence") or {}).get("value", "")
            protein = (((res.get("proteinDescription") or {}).get("recommendedName") or {})
                       .get("fullName") or {}).get("value")
            if not protein:
                subs = (res.get("proteinDescription") or {}).get("submissionNames") or []
                protein = (subs[0].get("fullName") or {}).get("value") if subs else None
            yield EstherHit(
                accession=res.get("primaryAccession", ""),
                protein_name=protein,
                organism=org.get("scientificName"),
                taxid=org.get("taxonId"),
                lineage="; ".join(lineage) if isinstance(lineage, list) else lineage,
                sequence=seq, length=len(seq), family=_family_of(res),
            )
            seen += 1
            if max_results is not None and seen >= max_results:
                return
        # The next page URL already carries the cursor: re-sending params restarts the
        # walk at page one, forever.
        m = _LINK_NEXT.search(resp.headers.get("Link", ""))
        url, params = (m.group(1) if m else None), None


def build_query(length_min: int, length_max: int, bacteria_only: bool = True,
                require_signal_peptide: bool = True) -> str:
    """Length-, taxonomy- and secretion-restricted ESTHER slice.

    The length window is derived from the positives rather than hardcoded, so negatives
    cannot be separable by length alone (spec section 5.2, first matching criterion).

    require_signal_peptide is a fourth matching axis the spec does not list but the data
    demanded. Without it the trivial baseline scored AUC 0.954 cluster-grouped, and the
    coefficients showed why: negatives were Leu-rich, positives Ser/Thr/Gly/Pro-rich,
    which is a secreted-versus-cytoplasmic signature. Every characterised polyesterase is
    secreted (IsPETase 1-27, LCC 1-34, TfCut2 1-40 all carry signal peptides), so a
    negative set drawn from mostly cytoplasmic families lets the model win by learning
    "is this exported?" instead of "does this degrade PET?".
    """
    parts = ["database:esther", f"length:[{length_min} TO {length_max}]"]
    if bacteria_only:
        parts.append("taxonomy_id:2")
    if require_signal_peptide:
        parts.append("(ft_signal:*)")
    return " AND ".join(parts)


def harvest(length_min: int, length_max: int, max_scan: int = 20000,
            bacteria_only: bool = True, require_signal_peptide: bool = True) -> Dict[str, object]:
    """Scan the ESTHER slice and bucket entries into negatives, near misses and rejects.

    Returns a report with the hits plus the family census, which is worth keeping: it
    shows what was passed over and why, and makes the allowlist auditable rather than
    something to take on faith.
    """
    negatives: List[EstherHit] = []
    near_misses: List[EstherHit] = []
    census: Counter = Counter()
    rejected: Counter = Counter()
    n_scanned = 0

    query = build_query(length_min, length_max, bacteria_only, require_signal_peptide)
    for hit in stream(query, max_results=max_scan):
        n_scanned += 1
        census[hit.family or "(no family)"] += 1

        if not hit.is_clean:
            rejected["non-standard residues"] += 1
            continue
        if hit.family in POSITIVE_FAMILIES:
            rejected["positive family"] += 1
            continue
        if hit.family in POLYESTER_ACTIVE_EXCLUDE:
            rejected["polyester-active family"] += 1
            continue

        role = hit.role
        if role == "negative":
            negatives.append(hit)
        elif role == "near_miss":
            near_misses.append(hit)
        else:
            rejected["family not on allowlist"] += 1

    return {
        "negatives": negatives,
        "near_misses": near_misses,
        "census": census,
        "rejected": rejected,
        "n_scanned": n_scanned,
    }
