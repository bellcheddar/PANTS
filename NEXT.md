# PANTS — where this is, and what to do next

Written so the project can be picked up cold. Live at
[pants.mdeller.com](https://pants.mdeller.com). Data and artefacts:
[10.5281/zenodo.21807353](https://doi.org/10.5281/zenodo.21807353) (concept DOI, always
resolves to the newest version; current is v0.3.0).

## What it is

A search for PET-degrading enzymes that work at **37 °C in serum** — a therapeutic framing —
rather than at 70 °C in an industrial reactor. It mines metagenomes, folds and measures
candidates, and ranks them against a reference set of characterised enzymes.

## What it found

Read `FINDINGS.md`. It is **generated** from JSON artefacts by
`scripts/build_findings_document.py`, so its numbers cannot go stale; do not edit it by hand.
In short:

1. Active-site geometry does not predict PET activity on measured labels (AUC 0.398, having
   been 0.808 at p 1.7×10⁻⁷ on labels inferred from database absence).
2. The determinants are **lineage-specific** — the same model fitted inside different
   families points in unrelated or opposite directions — which is why nothing transfers.
3. A learned classifier does not beat nearest-neighbour retrieval, so the shipped ranking
   **is** retrieval, with a stated competence band on every candidate.

The project's real contribution is that third point plus the negative results, not the
classifier it set out to build. Do not quietly re-add a classifier without new data.

## The one thing that would change any of this

`release/validation_panel.csv` — **150 assays across 15 new 30%-identity lineages**, one
protocol, 37 °C on amorphous PET. 80 are PANTS' own metagenome candidates (68 from human
gut, the therapeutically relevant environment) and 70 are homologues from the Science 2025
search space. It is designed to take **evaluable lineages from 3 to about 18**, which is the
binding constraint on every result above.

Selection is by breadth, deliberately not by predicted activity: using the prediction to
choose what to measure would select enzymes resembling what is already known, which is the
bias the project exists to correct. Regenerate with
`scripts/design_validation_panel.py`.

Two cheaper changes that also help: report quantitative values with explicit limits of
detection for negatives rather than binary calls, and cross the assay conditions so
within-enzyme temperature contrasts are identifiable.

## How to run it

```bash
# two virtual environments, deliberately
.venv/bin/python       # science: torch, biotite, sklearn
.venv-web/bin/python   # web: Flask + gunicorn ONLY -- never let torch in here

.venv/bin/python -m pytest tests/ -q      # 68 tests, all should pass
./deploy.sh                                # code + database + structures
./deploy.sh --no-data                      # code only, when the DB is being written
```

`deploy.sh` validates everything **before** the first rsync, including that the thin venv can
import the app. That ordering is deliberate: a failure after the code rsync leaves new code
on the droplet with the service running the old build and nothing saying so.

## Things that will bite you

- **Coordinate URLs must carry `?v=`.** nginx serves `/static/` as
  `max-age=31536000, immutable`, so a browser never revalidates. An unversioned structure URL
  meant users saw a cached multimer for weeks after the server was serving monomers. Use
  `coord_url()`; never build a `/static/structures/...` path by hand.
- **`request.args.get()` returns `None`, and `None` is a real key** in the lineage dict (the
  enzymes with no assigned wild type). Distinguish absence from value.
- **Bulk identifiers all contain a colon**; curated names never do. That is the rule the
  named-enzyme filter uses, because a denylist of prefixes was silently broken by every new
  bulk source.
- **`CREATE TABLE IF NOT EXISTS` never reaches an existing database.** New columns go in
  `COLUMN_MIGRATIONS`, new tables in `LATE_TABLES`.
- **Never pool predicted and experimental coordinates** in a geometric analysis. Paired within
  the same protein, the oxyanion hole differs at p 1.1×10⁻⁶. The `source` column exists for
  this.
- **Evaluate by cluster, never by sequence**, and prefer leave-one-cluster-out with
  size weighting. One cluster holds 68% of the labelled data, so k-fold splits are nearly
  forced and their seed-to-seed variance is an artefact — a `± 0.000` across ten seeds once
  meant two distinct partitions, not stability.
- **ESMFold is ~7 min per structure on this machine** and is the only expensive operation.
  Everything else is minutes.

## Deliberately not done

- **MHETase pipeline.** The second enzyme of the PET pathway, in the Tannase family, so a
  PETase-seeded profile search cannot reach it. Needs its own seeds, profiles and negatives.
  It is scope expansion, not a fix for anything above.
- **PU learning on the 25,041 unlabelled homologues.** Its assumptions are violated: the
  positives were selected by researchers who study thermophilic cutinases, not at random, so
  it would import collection bias as signal and trade away the 150 measured negatives — the
  rarest asset here — for pseudo-labels.
- **Domain-adversarial training.** No model-selection criterion is available on three
  evaluable folds, and since lineage partly predicts label, a successful adversary destroys
  signal and confound together with no way to tell which happened.
