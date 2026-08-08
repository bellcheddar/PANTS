"""Guard the findings against silent drift.

The three results this project reports are read from JSON artefacts at build time rather
than typed into prose, which stops them going stale. It does not stop something subtler: an
artefact being regenerated with a changed shape, or an analysis quietly starting to select
different enzymes, so the document renders cleanly with numbers that no longer mean what the
surrounding sentences claim.

These tests assert the STRUCTURE and the DIRECTION of each finding, never the exact value.
A test pinning `grouped_auc_mean == 0.398` would fail the first time anyone legitimately
adds data, which trains people to delete tests. A test asserting that the measured-label
geometry result stays below the inferred-label one it superseded fails only when the finding
itself has changed, which is the moment somebody should look.
"""

from __future__ import annotations

import json
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
RELEASE = ROOT / "release"
INTERIM = ROOT / "data" / "interim"


def artefact(name: str):
    for base in (RELEASE, INTERIM):
        p = base / name
        if p.exists():
            return json.loads(p.read_text())
    pytest.skip(f"{name} not built yet")


def test_geometry_on_measured_labels_is_at_or_below_chance():
    """Finding 1. Geometry gave 0.808 on inferred labels; on measured ones it must not."""
    g = artefact("geometry_measured_labels.json")
    assert g["n_positive"] > 100 and g["n_negative"] > 100
    assert g["coordinate_source"] == "esmfold only, both classes", (
        "mixing coordinate sources reintroduces the confound this analysis controls for")
    assert g["grouped_auc_mean"] < 0.55, (
        "geometry now separates the measured classes; the reported finding has changed")
    assert g["per_feature"]["cleft_depth_A"]["p"] > 0.01, (
        "cleft depth has become significant on measured labels; re-examine Finding 1")


def test_lineages_disagree_more_than_bootstrap_noise():
    """Finding 2. Between-lineage agreement must stay well below the within-lineage ceiling."""
    lin = artefact("p0_lineage_determinants.json")
    for features, r in lin.items():
        ceiling = sum(r["ceiling"].values()) / len(r["ceiling"])
        between = [v["mean"] for v in r["between"].values()]
        assert ceiling > 0.4, (
            f"{features}: the within-lineage ceiling has collapsed, so the test no longer "
            f"has the power to distinguish disagreement from noise")
        assert max(between) < ceiling - 0.2, (
            f"{features}: lineages now agree about as well as replicates of one lineage; "
            f"the lineage-specificity finding no longer holds")


def test_bigger_language_models_do_not_help():
    """Finding 3, part one. Capacity was never the bottleneck."""
    v = artefact("sequence_head_variants.json")
    sizes = [k for k in v if k.startswith("esm2-")]
    assert len(sizes) >= 2, "need at least two model sizes to make the comparison"
    aucs = {k: v[k]["auc_size_weighted"] for k in sizes if "auc_size_weighted" in v[k]}
    assert max(aucs.values()) - min(aucs.values()) < 0.25, (
        "model size now changes the result substantially; Finding 3 needs revisiting")


def test_every_candidate_carries_a_competence_band():
    """Finding 3, part two. A ranked list without a stated competence range is the failure
    mode this project argues against, so the band is not optional."""
    r = artefact("candidate_retrieval_scores.json")
    counts = r["counts"]
    assert set(counts) >= {"in-range", "marginal", "out-of-range"}
    assert sum(counts.values()) == len(r["rows"])
    assert counts["out-of-range"] > 0, (
        "no candidate is out of range, which would mean the competence band stopped being "
        "computed rather than that the catalogue improved")
    for row in r["rows"]:
        assert row["band"] in {"in-range", "marginal", "out-of-range", "no match"}


def test_findings_document_has_no_unresolved_values():
    """The document is generated; a missing artefact must fail loudly, not render 'nan'."""
    doc = ROOT / "FINDINGS.md"
    if not doc.exists():
        pytest.skip("FINDINGS.md not built yet")
    import re

    text = doc.read_text()
    # Bare-substring matching is wrong here and this test proved it by failing on the "nan"
    # inside "determinants". Match standalone tokens only.
    for pattern, what in ((r"(?<![A-Za-z])nan(?![A-Za-z])", "nan"),
                          (r"(?<![A-Za-z])None(?![A-Za-z])", "None"),
                          (r"[-+]?inf(?![A-Za-z])", "inf")):
        hit = re.search(pattern, text)
        assert hit is None, (
            f"FINDINGS.md contains an unresolved {what} at offset {hit.start()}; "
            f"an artefact is missing or changed shape")
    # An unrendered f-string placeholder, as distinct from a markdown code span.
    assert not re.search(r"\{[a-z_]+[\[\.']", text), (
        "FINDINGS.md contains an unrendered template placeholder")
