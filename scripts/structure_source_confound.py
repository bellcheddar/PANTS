#!/usr/bin/env python3
"""Is active-site geometry comparable between a crystal structure and a model?

This started as a puzzle in the numbers. Re-running the geometry-versus-activity test on
the finished set gave a cluster-grouped AUC of 0.553 across 342 enzymes, but 0.749 on the
287 predicted ones alone -- adding 55 crystal structures made the result WORSE, stably,
with non-overlapping ranges across twenty random split seeds. More data lowering a score
is a confound announcing itself.

Among enzymes of the same activity class, geometry alone tells a crystal structure from a
model at cluster-grouped AUC 0.723. That could still be an artefact of which enzymes get
crystallised: the deposited ones are the famous, heavily engineered ones, and they might
genuinely differ. The only comparison that removes that is paired -- the same protein
measured both ways -- which is what this does: for every enzyme now carrying a deposit,
fetch its AlphaFold model as well, superpose and measure it through exactly the same
geometry.measure(), and compare feature by feature with a Wilcoxon signed-rank test.

The answer decides how the evaluation may be run. If the features are source-invariant,
experimental and predicted structures can be pooled and the 55 deposits are free extra
data. If they are not, pooling is a methodological error and any geometry model must
either hold the source constant or use only the features that survive.
"""

from __future__ import annotations

import json
import pathlib
import sys
from typing import Dict, List

import numpy as np
from scipy import stats

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config
from pipeline.db import connect
from pipeline.structure import fold, geometry, reference

CACHE = config.INTERIM_DIR / "paired_af"

# The reference_geometry column and the ActiveSite attribute holding the same quantity are
# NOT identically named ("_dist_A" against "_A"). Written out rather than derived, because
# a getattr on a name that does not exist returns None and silently drops the feature from
# the comparison instead of failing -- which it did, losing three of the nine on first run.
COLUMNS: List[tuple] = [
    ("ser_og_his_ne2_dist_A", "ser_og_his_ne2_A"),
    ("his_nd1_asp_od_dist_A", "his_nd1_asp_od_A"),
    ("ser_his_asp_angle_deg", "ser_his_asp_angle_deg"),
    ("oxyanion_n1_dist_A", "oxyanion_n1_A"),
    ("oxyanion_n2_dist_A", "oxyanion_n2_A"),
    ("oxyanion_n2_angle_deg", "oxyanion_n2_angle_deg"),
    ("cleft_width_A", "cleft_width_A"),
    ("cleft_depth_A", "cleft_depth_A"),
    ("n_cleft_residues", "n_cleft_residues"),
]


def _safe(enzyme_id: str) -> str:
    """The same filename rule reference.build uses, so the cache keys line up."""
    return enzyme_id.replace("/", "_").replace("*", "s")


def build_models() -> Dict[str, geometry.ActiveSite]:
    """AlphaFold model per deposit-carrying enzyme, superposed and measured identically."""
    CACHE.mkdir(parents=True, exist_ok=True)
    ref_cif = config.STATIC_DIR / "reference" / f"{config.ISPETASE_REFERENCE_PDB}.cif"
    with connect() as c:
        rows = c.execute(
            "SELECT r.enzyme_id, e.uniprot FROM reference_structures r "
            "JOIN characterised_enzymes e USING(enzyme_id) "
            "WHERE r.source='pdb' AND e.uniprot IS NOT NULL AND e.uniprot != ''").fetchall()

    sites: Dict[str, geometry.ActiveSite] = {}
    for enzyme_id, acc in rows:
        dest = CACHE / f"{_safe(enzyme_id)}.pdb"
        if not dest.exists():
            text = reference._fetch_alphafold(acc.split("-")[0])
            if not text:
                continue
            cif, _, _ = fold.superpose_onto_reference(
                reference.keep_first_polymer_chain(text), ref_cif)
            if cif is None or fold.write_viewer_pdb(cif, dest) is None:
                continue
        site = geometry.measure(dest)
        if site.ser_resnum is not None:
            sites[enzyme_id] = site
    return sites


def main() -> int:
    sites = build_models()
    with connect() as c:
        xray = {r["enzyme_id"]: dict(r) for r in c.execute(
            "SELECT enzyme_id, " + ", ".join(col for col, _ in COLUMNS)
            + " FROM reference_geometry")}

    print(f"paired on {len(set(sites) & set(xray))} enzymes "
          f"(same protein, crystal structure and AlphaFold model)\n")
    print(f"{'feature':24s} {'x-ray':>8s} {'model':>8s} {'delta':>8s}  {'paired p':>9s}      n")

    out = {}
    for col, attr in COLUMNS:
        a, b = [], []
        for enzyme_id, site in sites.items():
            v_x = (xray.get(enzyme_id) or {}).get(col)
            v_m = getattr(site, attr, None)
            if v_x is None or v_m is None:
                continue
            a.append(float(v_x))
            b.append(float(v_m))
        if len(a) < 8:
            print(f"{col:24s} too few pairs ({len(a)})")
            continue
        a, b = np.array(a), np.array(b)
        p = float(stats.wilcoxon(a, b).pvalue)
        star = "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else "ns"
        print(f"{col:24s} {a.mean():8.2f} {b.mean():8.2f} {a.mean() - b.mean():+8.2f}  "
              f"{p:9.2g} {star:4s} {len(a)}")
        out[col] = {"xray": a.mean(), "model": b.mean(), "delta": a.mean() - b.mean(),
                    "p": p, "n": len(a), "source_invariant": p >= 0.05}

    dest = config.INTERIM_DIR / "structure_source_confound.json"
    dest.write_text(json.dumps(out, indent=2, default=float))
    inv = [k for k, v in out.items() if v["source_invariant"]]
    print(f"\nsource-invariant at p>=0.05: {', '.join(inv) if inv else 'none'}")
    print(f"wrote {dest}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
