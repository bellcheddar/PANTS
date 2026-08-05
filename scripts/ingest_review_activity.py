"""Ingest the optimum temperatures and headline performance claims for the named variants.

These are the values the README's therapeutic-gap table already shows, moved out of prose
and into the database so the web app renders them with their provenance attached instead
of hardcoding numbers in a template.

**Provenance is deliberately weaker than the UniProt block, and marked as such.** These
come from a literature review rather than from reading each primary paper, so they carry
`ECO:0000305` (inferred by curator) and `extraction_confidence='review'`, against the
`ECO:0000269` (experimental, from a publication) on the UniProt-derived rows. Anything
reported from this set must be separable from the experimental set, which is why the code
is stored per row rather than assumed.

**Ranges keep their range.** Several are published as an interval (60 to 65 °C). The
numeric column takes the midpoint, because a table needs something sortable, and the
verbatim published wording goes in `raw_text` so the interval is never lost and can always
be shown instead of the midpoint. A midpoint presented as if it were a measured single
value would be a small fabrication.

IsPETase is deliberately given a SECOND row rather than having its existing one replaced.
UniProt curates 40 °C, this review gives 30 to 35, they were measured on different
substrates under different assays, and the README states both rather than picking a
winner. Overwriting would silently resolve a disagreement the project has chosen to show.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline.db import connect, retry_write

# The PRIMARY paper for each enzyme, not the review these values were collated from.
# A reader following a citation should land on the work that made the measurement; a
# secondary write-up is a route to the primary source, not a substitute for it. Every DOI
# here was verified through Crossref: resolved, and the type, title, journal and year
# checked against the citation.
LEGACY_SOURCE = "https://marcdeller.com/engineering-evolution-how-fast-petase-and-other-variants-are-transforming-plastic-biodegradation/"

# enzyme_id, midpoint °C, verbatim wording, headline performance claim, primary DOI
ROWS = [
    ("Z1-PETase",    30.0, "30 °C",              "13 mutations, two engineered disulfides; 40x expression yield", "10.1016/j.jhazmat.2023.132297"),
    ("IsPETase",     32.5, "30 to 35 °C",        "Wild type; weak on crystalline PET",                             "10.1126/science.aad6359"),
    ("DuraPETase",   37.0, "37 °C",              "10 mutations; +31 °C thermostability, ~300x activity",           "10.1021/acscatal.0c05126"),
    ("FAST-PETase",  50.0, "50 °C",              "38x activity; 33.8 mM monomers in 96 h",                         "10.1038/s41586-022-04599-z"),
    ("DepoPETase",   50.0, "~50 °C (applied)",   "7 mutations; melting temperature +23.3 °C, ~1407x product",      "10.1016/j.xcrp.2024.102295"),
    ("HotPETase",    62.5, "60 to 65 °C",        "21 mutations; melting temperature 82.5 °C",                      "10.1038/s41929-022-00821-3"),
    ("Cut190**SS",   65.0, "65 °C",              "Calcium-dependent conformational switching",                     "10.1021/acs.biochem.8b00624"),
    ("TurboPETase",  66.5, "65 to 68 °C",        "98.2% depolymerisation at 200 g/kg in 8 h",                      "10.1038/s41467-024-45662-9"),
    ("LCC-ICCG",     68.5, "65 to 72 °C",        "1.3 g PET waste in 3 days from 1.25 mg enzyme",                  "10.1038/s41586-020-2149-4"),
    ("LCC-A2",       78.0, "78 °C",              "LCC-ICCG plus H218Y/N248D",                                      "10.1002/pro.70282"),
    ("ThermoPETase", 50.0, "50 °C",              "3 mutations; the thermostabilised scaffold FAST-PETase was built on", "10.1021/acscatal.9b00568"),
]


def main() -> None:
    added = skipped = 0
    with connect() as c:
        known = {e for e in c.execute("SELECT enzyme_id FROM characterised_enzymes")
                 for e in [e[0]]}
        # Rows written by an earlier run pointed at the review rather than the primary
        # paper. Repoint them rather than inserting duplicates.
        for enzyme_id, _m, _v, _c, doi in ROWS:
            c.execute("UPDATE activity_measurements SET source_doi=? "
                      "WHERE enzyme_id=? AND source_doi=?",
                      (doi, enzyme_id, LEGACY_SOURCE))
        existing = {(r[0], r[1]) for r in c.execute(
            "SELECT enzyme_id, parameter_type FROM activity_measurements "
            "WHERE extraction_confidence = 'review'")}

    for enzyme_id, midpoint, verbatim, claim, doi in ROWS:
        if enzyme_id not in known:
            print(f"  SKIP {enzyme_id}: not in characterised_enzymes")
            skipped += 1
            continue

        def _do(enzyme_id=enzyme_id, midpoint=midpoint, verbatim=verbatim, claim=claim, doi=doi):
            with connect() as c:
                if (enzyme_id, "topt") not in existing:
                    c.execute(
                        "INSERT INTO activity_measurements (enzyme_id, parameter_type, "
                        " substrate_form, temperature_c, rate_value, rate_units, raw_text, "
                        " evidence_code, source_doi, extraction_confidence) "
                        "VALUES (?,?,?,?,?,?,?,?,?,?)",
                        (enzyme_id, "topt", "PET", midpoint, midpoint, "degC",
                         f"Published as {verbatim}."
                         + (" Midpoint stored in the numeric column; the published "
                            "interval is this text." if " to " in verbatim else ""),
                         "ECO:0000305", doi, "review"))
                if (enzyme_id, "performance_claim") not in existing:
                    c.execute(
                        "INSERT INTO activity_measurements (enzyme_id, parameter_type, "
                        " substrate_form, raw_text, evidence_code, source_doi, "
                        " extraction_confidence) VALUES (?,?,?,?,?,?,?)",
                        (enzyme_id, "performance_claim", "PET", claim,
                         "ECO:0000305", doi, "review"))
        retry_write(_do)
        added += 1

    with connect() as c:
        n = c.execute("SELECT COUNT(*) FROM activity_measurements "
                      "WHERE extraction_confidence='review'").fetchone()[0]
        blog = c.execute("SELECT COUNT(*) FROM activity_measurements "
                         "WHERE source_doi LIKE '%marcdeller.com%'").fetchone()[0]
        print(f"  rows still citing the blog: {blog}")
        tot = c.execute("SELECT COUNT(*) FROM activity_measurements").fetchone()[0]
    print(f"\n  processed {added}, skipped {skipped}")
    print(f"  rows from this review: {n} of {tot} measurements total")


if __name__ == "__main__":
    main()
