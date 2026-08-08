# `scripts/` — what each one does

Run everything with `.venv/bin/python` (the science environment). `.venv-web` carries Flask
and gunicorn only and will not import these.

Ordered by what you are trying to do, not alphabetically.

## Understand what this project found

| Script | Does |
|---|---|
| `build_findings_document.py` | Regenerates `FINDINGS.md`, reading every figure from the JSON artefact that produced it. **`FINDINGS.md` is generated — edit this, not that.** |
| `p0_lineage_specific_determinants.py` | The go/no-go test. Fits the same model inside each lineage and compares directions against a bootstrap ceiling. Finding 2. |
| `geometry_measured_labels.py` | Geometry versus activity on measured labels. Finding 1. |
| `identity_decay_curve.py` | Where the signal holds, as a function of sequence identity, with paired bootstrap intervals. |
| `sequence_head_variants.py` | ESM-2 35M/150M/650M compared under one leave-one-cluster-out evaluation. |
| `structure_source_confound.py` | Paired crystal-versus-model geometry on the same protein. The control behind "never pool coordinate sources". |
| `run_evaluation_protocol.py` | The full evaluation protocol, including the components that are not evaluable and why. |
| `eval_within_family.py` | The within-family contrast across every negative definition. |
| `geometry_vs_activity.py` | Superseded by `geometry_measured_labels.py`; kept because it produced the 0.808 result the newer one overturned. |
| `geometry_vs_activity_v2.py` | The intermediate re-run on the finished structure set, before measured labels existed. |
| `active_site_embeddings.py` | Pools the language model over active-site residues only. **Written, never run** — it would be measured on the same ten lineages, so it could not produce a trustworthy answer. |

## Decide what to do next

| Script | Does |
|---|---|
| `design_validation_panel.py` | The 150-assay panel: 15 new lineages, selected by breadth rather than by prediction. Writes `release/validation_panel.csv`. |
| `score_candidates_by_retrieval.py` | The shipped ranking. Identity to the nearest **measured-active** enzyme, with a competence band per candidate. |

## Build and maintain the data

| Script | Does |
|---|---|
| `run_recall.py` | Profile search over metagenome proteins, then the triad filter. |
| `run_embed.py` | ESM-2 embeddings for the training set and candidates. |
| `run_train.py` | The original three-scheme head. Superseded by `run_evaluation_protocol.py`. |
| `fold_reference_bulk.py` | ESMFold/AlphaFold for reference enzymes lacking a structure. ~7 min each; resumable. |
| `run_fold.py` | Folds candidates rather than reference enzymes. |
| `fold_drain.py` | Works through the deferred fold queue. |
| `rebuild_deposit_linked.py` | Rebuilds structures for enzymes that just acquired a PDB deposit. |
| `extract_primary_dois.py` | Lifts the citation out of free text into `primary_doi`. |
| `fetch_entry_dates.py` | UniProt first-public dates, for the prospective holdout. |
| `ingest_science_supplements.py` | Science 2025 S1 (homologues) and S7 (cluster ecology); documents why S2 is skipped. |
| `paper_access_audit.py` | Which multi-enzyme papers are openly reachable, and which need institutional access. |
| `promote_named_enzyme.py` | Renames a bulk-imported enzyme to its published name, cascading through the tables that reference it. |
| `repair_screen_ingest.py` | One-off repair of provenance and contradictory labels from the first screen ingest. Kept as the record of what was fixed. |
| `compute_lineage_identity.py` | Identity of each variant to its lineage wild type. |
| `curate_ec.py`, `nearmiss.py`, `retest.py`, `run_gut2.py`, `more_gut.py` | Reference-set curation and harvesting. |
| `reconcile_topt.py`, `set_headlines.py`, `harmonise.py` | Reference-set tidying for display. |
| `mutation_triad_distance.py`, `remeasure_geometry.py` | Active-site measurements on the reference structures. |

## Release and deploy

| Script | Does |
|---|---|
| `build_release.py` | Exports the CSVs, **recomputes every published figure from them**, and generates `DATASHEET.md`. Nothing published is transcribed. |
| `ingest_review_activity.py` | Curated measurements from the review literature. |

## Conventions that matter

- Every script inserts the repo root on `sys.path` — Python puts the *script's* directory
  there, not the working directory, so `import pipeline` fails without it.
- Anything that writes to the database opens a `stage_manifest`, so a run records its inputs,
  outputs, tool versions, git commit and wall time even if it later raises.
- Analyses write a JSON artefact next to their output. The site and `FINDINGS.md` read those
  artefacts rather than restating numbers, which is how a superseded AUC once survived in
  prose for weeks.
