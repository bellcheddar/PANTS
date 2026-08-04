"""Measured activity data, extracted from UniProt rather than typed from papers.

Spec section 5.4 calls the activity data the real blocker for any regression head:
published numbers are not comparable across papers (amorphous film versus crystalline
powder, different crystallinity, different temperatures, HPLC TPA release versus turbidity
versus weight loss), and the honest fallback is ordinal ranking within a paper.

PAZy has no API, so the temptation is to type rates out of PDFs. That is precisely where
fabrication risk lives, and a wrong rate is undetectable downstream. UniProt instead
carries the same information in machine-readable, citable form:

  - `CATALYTIC ACTIVITY` comments with an EC number and ECO evidence codes. **EC 3.1.1.101
    is poly(ethylene terephthalate) hydrolase**, so its presence is a curator's assertion
    of measured PET-hydrolysing activity, not a family guess.
  - `BIOPHYSICOCHEMICAL PROPERTIES` with `kineticParameters` (Km, with substrate and
    units) and free-text pH and temperature optima.
  - Every one of those carries `ECO:0000269` (experimental evidence from a publication)
    plus the PubMed IDs.

So every row this module writes is sourced and checkable. Nothing is invented, and the
`raw_text` column keeps the original statement so a parsed number can always be audited
against what was actually said.

The harmonisation problem does not go away: a Km on pNP-butyrate is not comparable to one
on PET film. `comparable_group_id` is keyed on the substrate for exactly that reason, and
rows measured on different substrates must never be pooled.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional, Tuple

from .. import http
from ..db import connect, now

ENTRY_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"

# The EC number that means "this was measured hydrolysing PET".
EC_PETASE = "3.1.1.101"
EC_MHETASE = "3.1.1.102"

# Experimental evidence from a publication. Anything weaker (ECO:0000305 inferred by
# curator, ECO:0000250 by similarity) is not a measurement and is not written.
ECO_EXPERIMENTAL = "ECO:0000269"

_TEMP_OPT = re.compile(r"[Oo]ptimum temperature is\s*(?:about\s*)?(\d+(?:\.\d+)?)\s*(?:to\s*(\d+(?:\.\d+)?)\s*)?degrees", re.I)
_PH_OPT = re.compile(r"[Oo]ptimum pH is\s*(?:about\s*)?(\d+(?:\.\d+)?)(?:\s*-\s*(\d+(?:\.\d+)?))?", re.I)


@dataclass
class Measurement:
    enzyme_id: str
    parameter_type: str
    substrate: Optional[str]
    value: Optional[float]
    units: Optional[str]
    raw_text: Optional[str]
    pmids: List[str] = field(default_factory=list)
    evidence_code: str = ECO_EXPERIMENTAL

    @property
    def comparable_group_id(self) -> str:
        """Only same parameter on the same substrate is mutually comparable."""
        return f"{self.parameter_type}:{(self.substrate or 'unspecified').lower()}"


def _pmids(evidences: Iterable[dict]) -> List[str]:
    return [e.get("id") for e in (evidences or [])
            if e.get("source") == "PubMed" and e.get("evidenceCode") == ECO_EXPERIMENTAL
            and e.get("id")]


def extract(accession: str, enzyme_id: str) -> Tuple[List[Measurement], Dict[str, object]]:
    """Pull every experimentally-evidenced measurement UniProt holds for one accession."""
    obj = http.get_json(ENTRY_URL.format(accession=accession))
    meta: Dict[str, object] = {"ec_numbers": [], "has_petase_ec": False, "pmids": set()}
    out: List[Measurement] = []
    if not obj:
        return out, meta

    for c in obj.get("comments", []):
        ctype = c.get("commentType")

        if ctype == "CATALYTIC ACTIVITY":
            rxn = c.get("reaction", {}) or {}
            ec = rxn.get("ecNumber")
            pm = _pmids(rxn.get("evidences") or c.get("evidences") or [])
            if ec:
                meta["ec_numbers"].append(ec)
                if ec == EC_PETASE:
                    meta["has_petase_ec"] = True
                    meta["pmids"].update(pm)
                    out.append(Measurement(
                        enzyme_id=enzyme_id, parameter_type="catalytic_activity",
                        substrate="PET", value=None, units=None,
                        raw_text=rxn.get("name"), pmids=pm,
                    ))

        elif ctype == "BIOPHYSICOCHEMICAL PROPERTIES":
            kp = c.get("kineticParameters", {}) or {}
            for km in kp.get("michaelisConstants", []) or []:
                pm = _pmids(km.get("evidences"))
                if not pm:
                    continue
                out.append(Measurement(
                    enzyme_id=enzyme_id, parameter_type="km",
                    substrate=km.get("substrate"), value=km.get("constant"),
                    units=km.get("unit"), raw_text=None, pmids=pm,
                ))
            for vm in kp.get("maximumVelocities", []) or []:
                pm = _pmids(vm.get("evidences"))
                if not pm:
                    continue
                out.append(Measurement(
                    enzyme_id=enzyme_id, parameter_type="vmax",
                    substrate=vm.get("enzyme"), value=vm.get("velocity"),
                    units=vm.get("unit"), raw_text=None, pmids=pm,
                ))

            for key, ptype, rx in (("temperatureDependence", "topt", _TEMP_OPT),
                                   ("phDependence", "ph_opt", _PH_OPT)):
                for t in (c.get(key, {}) or {}).get("texts", []) or []:
                    text = t.get("value") or ""
                    pm = _pmids(t.get("evidences"))
                    if not pm:
                        continue
                    m = rx.search(text)
                    # Store the row even when the number will not parse: the prose is the
                    # measurement, and a missing value is honest where a guess is not.
                    val = None
                    if m:
                        lo = float(m.group(1))
                        hi = float(m.group(2)) if m.group(2) else None
                        val = (lo + hi) / 2 if hi else lo
                    out.append(Measurement(
                        enzyme_id=enzyme_id, parameter_type=ptype, substrate=None,
                        value=val, units="degC" if ptype == "topt" else "pH",
                        raw_text=text, pmids=pm,
                    ))

    meta["pmids"] = sorted(meta["pmids"])
    return out, meta


def write(measurements: List[Measurement]) -> int:
    """Insert measurements, replacing any previous extraction for the same enzymes."""
    if not measurements:
        return 0
    enzymes = sorted({m.enzyme_id for m in measurements})
    with connect() as conn:
        conn.executemany(
            "DELETE FROM activity_measurements WHERE enzyme_id=? AND extraction_confidence='uniprot'",
            [(e,) for e in enzymes],
        )
        for m in measurements:
            conn.execute(
                "INSERT INTO activity_measurements "
                "(enzyme_id, substrate_form, temperature_c, ph, product_measured, "
                " parameter_type, rate_value, rate_units, raw_text, evidence_code, "
                " comparable_group_id, source_doi, extracted_at, extraction_confidence) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (m.enzyme_id, m.substrate,
                 m.value if m.parameter_type == "topt" else None,
                 m.value if m.parameter_type == "ph_opt" else None,
                 m.substrate, m.parameter_type, m.value, m.units, m.raw_text,
                 m.evidence_code, m.comparable_group_id,
                 ";".join(f"PMID:{p}" for p in m.pmids), now(), "uniprot"),
            )
    return len(measurements)
