"""Choose ONE optimum temperature per enzyme, by a stated rule, and record why.

UniProt often curates several optima for the same enzyme from different papers and
different substrates, all in one free-text block. Something has to pick the number the
table shows, and until now that was positional accident: whichever value happened to be
parsed first. That produced two real errors.

**The rule, in order:**

1. Prefer a value measured on **PET** over one measured on a model ester. This project is
   about PET degradation, and pNP-butanoate is not PET. LCC was the casualty: 50 °C is its
   optimum on pNP-butanoate, while UniProt's own text says the optimum on PET is "superior
   to 70 degrees Celsius". Shown at 50 °C it sat near the therapeutic end of the table,
   which is the opposite of true for the substrate the site is about.
2. Failing that, prefer the value with the most independent publications behind it. TfCut2
   had four candidates; 65 °C carries three PMIDs and the stored 55 °C carried one.
3. Record the substrate and the alternatives in raw_text, so the discarded values stay
   visible rather than being silently resolved away.

A value that is a BOUND rather than an optimum keeps that distinction in its text, because
"superior to 70" is not the same claim as "70".
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline.db import connect, retry_write

# enzyme_id -> (value, replacement raw_text)
FIXES = {
    "LCC": (70.0,
        "On PET the optimum is SUPERIOR TO 70 degrees Celsius (PubMed:22194294); 70 is "
        "stored as a lower bound, not a measured optimum. The previously stored 50 degrees "
        "is the optimum on pNP-butanoate (PubMed:22194294, PubMed:24593046), a model ester "
        "rather than PET. Rule 1: prefer the PET substrate."),
    "TfCut2": (65.0,
        "Optimum temperature is 65 degrees Celsius (PubMed:15638529, PubMed:31690819, "
        "PubMed:32269349). UniProt also records 25-50 (PubMed:24728714), 55 "
        "(PubMed:23604968) and 60 (PubMed:20816933); none states a substrate. Rule 2: 65 "
        "carries three independent publications, the previously stored 55 carried one."),
    "PLC:G8GER6": (65.0,
        "Optimum temperature is 65 degrees Celsius (PubMed:32269349); UniProt also records "
        "55 (PubMed:23604968). Rule 2, applied for consistency with TfCut2."),
}

# TurboPETase came from the review pass, not UniProt, and its stated interval is unsupported
TURBO = (65.0,
    "The primary paper (doi:10.1038/s41467-024-45662-9) runs its depolymerisation at 65 "
    "degrees Celsius and never reports an optimum for TurboPETase: it tests 50, 60 and 65 "
    "and 65 performs best of those three. The previously stored '65 to 68' interval does "
    "not appear in it, and may have been conflated with the melting temperature of 84 "
    "degrees. Stored as the operating temperature, not a measured optimum.")


def main() -> None:
    def _do() -> None:
        with connect() as c:
            for eid, (val, text) in FIXES.items():
                c.execute(
                    "UPDATE activity_measurements SET rate_value=?, temperature_c=?, "
                    "raw_text=? WHERE enzyme_id=? AND parameter_type='topt' "
                    "AND evidence_code='ECO:0000269'", (val, val, text, eid))
            c.execute(
                "UPDATE activity_measurements SET rate_value=?, temperature_c=?, raw_text=? "
                "WHERE enzyme_id='TurboPETase' AND parameter_type='topt'",
                (TURBO[0], TURBO[0], TURBO[1]))
    retry_write(_do)

    with connect() as c:
        print(f"  {'enzyme':<14}{'Topt':>6}  basis")
        for eid in list(FIXES) + ["TurboPETase", "IsPETase"]:
            r = c.execute("SELECT rate_value, raw_text FROM activity_measurements "
                          "WHERE enzyme_id=? AND parameter_type='topt' "
                          "ORDER BY CASE evidence_code WHEN 'ECO:0000269' THEN 0 ELSE 1 END "
                          "LIMIT 1", (eid,)).fetchone()
            if r:
                print(f"  {eid:<14}{r[0]:>6}  {(r[1] or '')[:74]}")


if __name__ == "__main__":
    main()
