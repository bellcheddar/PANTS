"""One-line descriptions for the summary table.

Each is a condensation of the curated note this project already holds for that enzyme, not
a new claim. Enzymes with a published performance figure get that instead, from
activity_measurements, and are absent here.

These are descriptions, not measurements, which is why they are a column on
characterised_enzymes rather than rows in activity_measurements: "hydrolyses MHET to TPA
and EG" is what an enzyme IS, and storing it beside measured optimum temperatures with an
evidence code would misrepresent it as a result.
"""

from __future__ import annotations

import pathlib
import sys

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline.db import connect, retry_write

HEADLINES = {
    "IsPETase": "Wild type from Piscinibacter sakaiensis; weak on crystalline PET",
    "LCC": "Leaf-branch compost metagenome cutinase; parent of the industrial ICCG variant",
    "TfCut2": "Thermobifida fusca cutinase 2; thermophilic industrial lineage",
    "Cut190": "Saccharomonospora viridis AHK190 cutinase; calcium-dependent stabilisation",
    "BhrPETase": "Bacterium HR29 hydrolase; the parent TurboPETase was engineered from",
    "MHETase": "Second enzyme of the pathway, hydrolysing MHET to TPA and EG. Tannase group, not the PETase family",
    "IsPETase-W159H/S238F": "The narrowed-cleft double mutant that showed PETase's cleft is wider than a cutinase's",
    "ThermoPETase": "Thermostabilised IsPETase; the scaffold FAST-PETase was built on",
    "HGMP01": "Human gut metagenome; hydrolyses PET nanoparticles and outperformed every other gut candidate",
    "HGMP02": "Human gut metagenome polyesterase, assayed on PET nanoparticles",
    "HGMP03": "Human gut metagenome polyesterase, assayed on PET nanoparticles",
    "HGMP04": "Human gut metagenome polyesterase, assayed on PET nanoparticles",
    "HGMP05": "Human gut metagenome polyesterase, assayed on PET nanoparticles",
    "FAST-PETase-N212A/K233C/S282C": "Disulfide-stabilised FAST-PETase; C233-C282 crosslinks the same pair Z1-PETase uses",
}

if __name__ == "__main__":
    with connect() as c:
        known = {r[0] for r in c.execute("SELECT enzyme_id FROM characterised_enzymes")}
    missing = [k for k in HEADLINES if k not in known]
    for k in missing:
        print(f"  SKIP {k}: not in characterised_enzymes")

    def _do():
        with connect() as c:
            for eid, text in HEADLINES.items():
                if eid in known:
                    c.execute("UPDATE characterised_enzymes SET headline=? WHERE enzyme_id=?",
                              (text, eid))
    retry_write(_do)

    with connect() as c:
        blank = [r[0] for r in c.execute("""
            SELECT ce.enzyme_id FROM characterised_enzymes ce
            WHERE ce.enzyme_id NOT LIKE '%:%' AND ce.is_positive=1
              AND (ce.headline IS NULL OR ce.headline='')
              AND NOT EXISTS (SELECT 1 FROM activity_measurements a
                              WHERE a.enzyme_id=ce.enzyme_id
                                AND a.parameter_type='performance_claim')""")]
    print(f"  set {len(HEADLINES) - len(missing)} headlines")
    print(f"  named enzymes still with NO headline and no performance claim: {blank or 'none'}")
