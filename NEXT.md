# PANTS — full context for resuming

Everything needed to pick this up cold. Live at [pants.mdeller.com](https://pants.mdeller.com).
Data and artefacts: [10.5281/zenodo.21807353](https://doi.org/10.5281/zenodo.21807353)
(concept DOI — always resolves to the newest version; currently **v0.3.0**).

Last worked on 2026-08-08. Working tree clean, 68 tests passing, nothing running.

---

## 1. What the project is

A search for PET-degrading enzymes that work at **37 °C in serum** — a therapeutic framing —
rather than at 70 °C in an industrial reactor, which is what the field optimises for. It mines
metagenomes for polyesterase folds, folds and measures the candidates, and ranks them against
a reference set of characterised enzymes.

## 2. Current state

| | |
|---|---|
| Catalogue | **1,348** characterised enzymes |
| — measured active | 389 |
| — **measured inactive** | **150** (expressed, assayed, no product released) |
| Activity measurements | 2,759 (2,757 with a DOI) |
| Metagenome candidates | 439, from 14.8M sequences scanned |
| Structures | 559 reference (60 crystal, 104 AlphaFold, 395 ESMFold) + 416 candidate |
| Unlabelled homologue pool | 25,041 |
| Schema | v17 |
| Data version | 0.3.0 |

Site pages: Home, Catalogue, Superpose, Lineages, Reference set, Methods, Stats.

## 3. What it found — read `FINDINGS.md`

`FINDINGS.md` is **generated** by `scripts/build_findings_document.py`, reading every figure
from the JSON artefact that produced it. **Do not edit it by hand** — regenerate it. This
exists because a Methods paragraph once quoted a superseded AUC for weeks after the run that
produced it had been redone.

Three findings, none of them what the project set out to prove:

**1. Active-site geometry does not predict PET activity.** On negatives inferred from database
absence, cleft depth separated the classes at AUC 0.808, p 1.7×10⁻⁷. On **measured** negatives:
**AUC 0.398**, cleft depth 0.529 at p 0.38. As negatives accumulated the result moved *toward*
chance (0.507 → 0.498 → 0.459 → 0.398), which is a confound diluting rather than an
underpowered test.

**2. The determinants are lineage-specific.** Fitting the same model inside each large lineage
and comparing directions: two fits of the *same* lineage agree at ~**+0.70**; two *different*
lineages agree at ~zero and in 2 of 3 pairs point in **opposite** directions. Same pattern in
two feature sets sharing nothing (3D coordinates, and a language model over sequence). This is
why nothing transfers, and why an 18× larger model changed nothing — **there may be no single
global rule to learn.**

**3. A learned classifier does not beat nearest-neighbour retrieval.** So the shipped ranking
*is* retrieval, with a competence band on every candidate: **16 in range, 39 marginal, 384 out
of range.** The 384 are precisely the novel enzymes the project set out to find, and nothing
here can score them — stated on the site rather than hidden.

**Do not quietly re-add a classifier without new data.** Three independent lines say the same
thing.

## 4. The one thing that would change any of this

**`release/validation_panel.csv`** — 150 assays across **15 new 30%-identity lineages**, one
protocol, **37 °C on amorphous PET**. Takes evaluable lineages from **3 to ~18**, which is the
constraint behind every result above.

- **80 PANTS candidates** (68 human gut, 10 marine plastisphere, 2 compost) — these are the
  out-of-range candidates nothing can currently score, so measuring them adds lineages *and*
  tests whether the mining found anything real.
- **70 Science-S1 homologues** for breadth.

Selected by breadth, **deliberately not by predicted activity** — using the prediction to
choose what to measure would select enzymes resembling what is already known, which is the bias
the project exists to correct. Regenerate with `scripts/design_validation_panel.py`.

Two cheaper improvements that also help: report quantitative values with explicit limits of
detection for negatives rather than binary calls, and cross the assay conditions so
within-enzyme temperature contrasts are identifiable.

## 5. Options from here, with an honest read

| Option | Assessment |
|---|---|
| **Run the panel** | The only thing that lifts the constraint. Not a coding task. |
| **Manuscript** | `FINDINGS.md` is manuscript-shaped but not a manuscript. Three negative results with deposited artefacts and a designed follow-up is a publishable unit, and negative results are scarce here precisely because databases record what worked. This is where accumulated value is least realised. |
| **MHETase pipeline** | The only open To Do item. Second enzyme of the pathway, Tannase family, so a PETase-seeded profile search cannot reach it — needs its own seeds, profiles, negatives. Well-defined, but it expands scope and would inherit the same lineage-confounding limit. |
| **Stop** | Defensible. Catalogue, live tool that states its own competence, three findings, citable deposit, handover docs. |

My read: **there is no more software worth writing until somebody runs assays.**

## 6. How to run it

```bash
.venv/bin/python        # science: torch, biotite, sklearn, openpyxl
.venv-web/bin/python    # web: Flask + gunicorn ONLY — never let torch in here

.venv/bin/python -m pytest tests/ -q   # 68 tests
./deploy.sh                             # code + database + structures
./deploy.sh --no-data                   # code only, when the DB is being written
```

`deploy.sh` validates everything **before** the first rsync, including that the thin venv can
import the app. A failure after the code rsync would leave new code on the droplet with the
service running the old build and nothing saying so.

Key scripts (38 in `scripts/`):

| Script | Does |
|---|---|
| `design_validation_panel.py` | the 150-assay panel |
| `build_findings_document.py` | regenerates `FINDINGS.md` |
| `p0_lineage_specific_determinants.py` | Finding 2, the go/no-go test |
| `geometry_measured_labels.py` | Finding 1 |
| `identity_decay_curve.py` | where the signal holds |
| `sequence_head_variants.py` | model-size comparison |
| `score_candidates_by_retrieval.py` | the shipped ranking |
| `run_evaluation_protocol.py` | full evaluation protocol |
| `structure_source_confound.py` | crystal-vs-model control |
| `build_release.py` | release bundle + datasheet, all figures recomputed |

## 7. Traps that have already cost time here

- **Coordinate URLs must carry `?v=`.** nginx serves `/static/` as
  `max-age=31536000, immutable`, so browsers never revalidate. An unversioned structure URL
  meant users saw a cached multimer for weeks while the server served monomers. Use
  `coord_url()`; never hand-build a `/static/structures/...` path.
- **`request.args.get()` returns `None`, and `None` is a real key** in the lineage dict (the
  enzymes with no assigned wild type). Distinguish absence from value.
- **Bulk identifiers all contain a colon**; curated names never do. That is the named-enzyme
  filter's rule, because a denylist of prefixes was silently broken by every new bulk source.
- **`CREATE TABLE IF NOT EXISTS` never reaches an existing database.** New columns go in
  `COLUMN_MIGRATIONS`, new tables in `LATE_TABLES`.
- **Never pool predicted and experimental coordinates** in a geometric analysis. Paired within
  the same protein, the oxyanion hole differs at p 1.1×10⁻⁶. The `source` column exists for this.
- **Evaluate by cluster, never by sequence.** Prefer leave-one-cluster-out with size weighting:
  one cluster holds 68% of the labelled data, so k-fold splits are nearly forced. A `± 0.000`
  across ten seeds once meant *two distinct partitions*, not stability. An unweighted mean of
  0.625 was seven clusters of 2–6 enzymes averaged against three real values of 0.443, 0.429,
  0.574 — the size-weighted answer was 0.472.
- **Ingests must APPEND to `activity_substrate_notes`, never overwrite.** A screen ingest once
  destroyed the provenance of 189 rows; recovered from the published Zenodo deposit.
- **ESMFold is ~7 min per structure** on this machine and is the only expensive operation.
- **Never type an identifier that lives in `credentials.env`.** An ORCID written from memory
  once put a stranger on this deposit; it was caught before publishing.

## 8. Deliberately not done, with reasons

- **PU learning on the 25k unlabelled homologues.** Assumptions violated — positives were
  selected by researchers who study thermophilic cutinases, not at random — so it would import
  collection bias as signal and trade away the 150 measured negatives, the rarest asset here,
  for pseudo-labels.
- **Domain-adversarial training.** No model-selection criterion on three evaluable folds, and
  since lineage partly predicts label, a successful adversary destroys signal and confound
  together with no way to tell which happened.
- **Masked-language-model adaptation on the homologue pool.** Sharpens the representation along
  the phylogenetic axis — exactly where the confound lives.
- **More ESMFold structures.** Geometry has not earned 7 min/structure.
- **Ordinal/ranking objectives.** Assessed and rejected on data: only 3 clusters hold two
  enzymes whose rates differ by more than replicate noise, against 10 for the binary objective,
  so ranking narrows the base rather than widening it.
