# PANTS: PETase ANnotation and Triage System

**Project brief for planning. Author: Marc C. Deller, D.Phil.**
**Target deployment: `pants.mdeller.com`**
**Status: pre-implementation. This document is input to a planning pass, not a finished plan.**

---

## 1. Objective

Discover, triage and engineer PETase-like and MHETase-like enzymes **for therapeutic use**, meaning degradation of PET microplastics under physiological conditions rather than in an industrial reactor.

The deliverable is an interactive web app at `pants.mdeller.com` that lets a user:

1. Browse a precomputed, ranked catalogue of candidate polyesterases mined from metagenomic and genomic sequence space.
2. Filter and sort candidates on **therapeutic** axes, not industrial ones.
3. Inspect any candidate's predicted structure, catalytic triad geometry and substrate cleft in 3D.
4. Submit a sequence of their own and have it scored against the same models.
5. Run in-silico engineering: score point mutants, propose stabilising or activity-enhancing substitutions, and compare variants side by side.
6. Interactively review all currently known PETase and MHETas structures, sequences and activity profiles including FAST-PETase, IsPETase, DuraPETase, HotPETase, LCC-ICCG, Z1-PETase, TurboPETase, Cut190**SS, HGMP01, and *ALL* others you can find.

### 1.1 The therapeutic constraint, stated up front

Essentially the entire published PETase field optimises for industrial conditions: above PET's glass transition (roughly 70 °C), often at alkaline pH, over hours to days, with enzyme cost and thermostability as the dominant concerns. LCC-ICCG, HotPETase, DuraPETase and the thermostable lineage generally are all pointed that way.

Therapeutic use inverts most of that:

| Axis | Industrial optimum | Therapeutic requirement |
|------|--------------------|-------------------------|
| Temperature | 65 to 75 °C | Peak activity at 37 °C |
| pH | 8 to 9 | 7.2 to 7.4 |
| Substrate | Amorphous, pre-treated, high surface area | Highly crystalline aged microplastic and nanoplastic |
| Medium | Buffer | Serum, plasma proteins, lipids, physiological ionic strength |
| Stability concern | Thermal denaturation over days | Serum protease resistance, aggregation, clearance |
| Immunogenicity | Irrelevant | Central |
| Product handling | Recovered and recycled | TPA and EG must be tolerable at achievable local concentrations |

A naive homology search seeded on characterised PETases ranks candidates **towards the industrial optimum**, because that is where the well-annotated, heavily-published enzymes sit. Correcting for that is a large part of what PANTS is for. Every ranking surface in the app should make the therapeutic axis explicit and default to it.

---

## 2. Why this is not a homology search

Detection is solved. Discrimination is not.

A profile HMM built from characterised PET hydrolases returns thousands of α/β-hydrolase fold members sharing the Ser-His-Asp triad and the oxyanion hole. Almost none have meaningful activity on crystalline PET. Sequence identity to IsPETase is a weak predictor of PET activity, so E-value rank is close to uninformative about the property of interest.

The learned model exists to do six things retrieval structurally cannot:

1. **Use negatives.** Retrieval has no way to express "same fold, same triad, functionally dead". Thousands of hard negatives are available from ESTHER: characterised lipases and carboxylesterases with no polyester activity. The discriminative signal is a scattered, non-contiguous set of positions (shallow wide cleft, aromatic clamp, the S214/I218 pairing giving IsPETase its mobile W185 where most cutinases have a fixed residue). A profile cannot weight a combination like that; an embedding head can.
2. **Invert the ranking for 37 °C**, by predicting activity and Topt as separate heads and ranking on the therapeutic composite.
3. **Continue past discovery into engineering.** You cannot BLAST a point mutant. ESM-2 masked-token log-odds give zero-shot variant effect scores on the same representation.
4. **Rank on several axes at once** via a multi-task head: activity, Topt, expressibility, aggregation propensity, protease susceptibility.
5. **Reach targets the seed does not resemble.** MHETase is not in the PETase family (tannase/feruloyl esterase group, with a lid domain), so a PETase-seeded profile search never reaches it. Embedding space crosses that gap.
6. **Return calibrated probabilities**, so a wet-lab screening budget can be set honestly: order the top 40, expect roughly N hits.

**Architecture that follows: retrieval is the recall stage, the model is the precision stage.** HMMER or MMseqs2 casts the net exhaustively and interpretably; the model ranks within it on axes retrieval is blind to. This also yields a free baseline: if the model cannot beat E-value rank on held-out characterised enzymes, that is a cheap and publishable negative result.

---

## 3. Architecture

Follow the **AlphaFraud pattern**, which is proven on this hardware:

- Batch pipeline on a systemd timer, writing to a single SQLite file.
- Always-on Flask plus gunicorn behind nginx, reading the same SQLite file.
- Single DigitalOcean droplet, roughly 3.8 GB RAM, certbot for TLS.
- No cluster, no external database, no queue.

### 3.1 Hard hardware constraints

- **All heavy compute is offline and precomputed.** Metagenome scanning, ESM-2 embedding of millions of sequences, and structure prediction never run in a web request.
- **On-demand scoring must fit alongside gunicorn in under about 2 GB.** ESM-2 **t12-35M** (roughly 150 MB fp32) is fine on CPU for single-sequence scoring. ESM-2 650M is not: reserve it for offline batch use only, if at all.
- Structure prediction (ESMFold or Boltz-2 via BoltzMaker) is **offline only**. Precompute structures for the top candidates and ship the mmCIF files as static assets.
- Target: on-demand single-sequence scoring returns in under 10 seconds on CPU.

### 3.2 Split of responsibilities

| Layer | Where it runs | Output |
|-------|---------------|--------|
| Recall (HMMER/MMseqs2 over metagenomes) | Offline, batch | Candidate FASTA |
| Embedding and multi-task scoring | Offline, batch | Scores in SQLite |
| Structure prediction and pocket geometry | Offline, batch | mmCIF assets plus geometry table |
| Catalogue browse, filter, plot | Flask, live | JSON to Plotly |
| 3D inspection | Client, live | Mol* over static mmCIF |
| User sequence submission | Flask, live, ESM-2 35M on CPU | Scores plus percentile against catalogue |
| Variant scoring | Flask, live, masked-token log-odds | Heatmap data |

---

## 4. Data sources

| Source | Use |
|--------|-----|
| **PAZy** | Characterised plastic-degrading enzymes. The positive set. Small (order 10^2) |
| **ESTHER** | α/β-hydrolase superfamily, family assignments, characterised non-polyesterases. The hard negative set |
| **PDB** | Experimental structures for IsPETase, LCC, TfCut2, Cut190, MHETase and engineered variants. Ground truth geometry |
| **UniProt / UniRef** | Reference sequence space, taxonomy, annotation quality flags |
| **MGnify, JGI IMG/M, OceanDNA, Tara Oceans** | Metagenomic assemblies for mining. Prioritise landfill, compost, marine plastisphere, wastewater, plastic-associated |
| **Meltome Atlas, FireProtDB** | Transfer learning for thermostability where PAZy data is too thin |
| **AlphaFold DB** | Precomputed structures where a UniProt match exists |
| **Literature (extracted)** | Quantitative activity data. See section 5.4 |

---

## 5. Model and data design

### 5.1 Recall stage

Profile HMMs from characterised PET hydrolases plus the relevant ESTHER families. MMseqs2 for speed across large metagenomic assemblies, HMMER for the final sensitive pass. Record E-value and profile identity for every candidate: these become the retrieval baseline the model must beat.

Filter for a complete catalytic triad and a recognisable oxyanion hole before scoring. Report how many candidates are discarded at this step.

### 5.2 Hard negative construction

This is the single most important dataset decision. Draw negatives from ESTHER families that share the fold and triad but are characterised as having no polyester activity: bacterial lipases, carboxylesterases, acetylcholinesterase-adjacent families, hormone-sensitive lipase family members. Match negatives to positives on:

- Sequence length distribution
- Overall identity to the nearest positive (so the model cannot win on identity alone)
- Taxonomic breadth

Explicitly include **near misses**: cutinases with confirmed activity on soluble esters but not on PET. These are the examples that define the boundary.

### 5.3 Positive-unlabelled learning

"Not annotated as a PETase" is not "tested and inactive". Most of sequence space is unlabelled, not negative. Treat the unlabelled pool as PU rather than negative: use a PU-aware loss (Elkan-Noto style class prior correction, or a bagging PU approach), and estimate the class prior explicitly rather than assuming it.

### 5.4 The activity data harmonisation problem

Published activity numbers are not comparable across papers: amorphous film versus crystalline powder versus nanoparticle, different crystallinity indices, different temperatures, different product quantitation (HPLC TPA release versus turbidity versus weight loss). This is a real blocker for any regression head.

Build a small extraction and harmonisation step: pull (enzyme, substrate form, crystallinity, temperature, pH, buffer, assay duration, product measured, rate) from the PAZy-linked literature into a table, and record which numbers are mutually comparable. Where harmonisation is impossible, **fall back to ordinal ranking within a paper** rather than pretending the absolute numbers are commensurate.

Flag clearly in the app which candidates are scored against harmonised quantitative data and which are scored against binary annotation only.

### 5.5 Model heads

Freeze ESM-2 t12-35M, train shallow heads. Do not fine-tune end to end: with order 10^2 positives, end-to-end fine-tuning will memorise.

| Head | Type | Training data |
|------|------|---------------|
| PET activity | Binary, calibrated | PAZy positives vs ESTHER hard negatives, PU-corrected |
| MHET activity | Binary, calibrated | MHETase family, separate seed |
| Topt | Regression | Meltome Atlas transfer, PAZy where available |
| Expressibility (soluble in *E. coli*) | Binary | TargetTrack, eSOL, or equivalent |
| Aggregation propensity | Regression | Published aggregation datasets, or TANGO/AGGRESCAN as weak labels |
| Therapeutic composite | Derived | Weighted combination, weights exposed in the UI |

The composite must be **transparent and user-adjustable in the app**: sliders for how much each axis matters, with the ranking updating live. This is more honest than a fixed black-box score and much more useful.

### 5.6 Engineering module

Zero-shot variant effect scoring from ESM-2 masked-token log-odds, plus:

- Per-position conservation and catalytic constraint, so the UI can grey out positions that should not be touched.
- Known beneficial substitutions from the literature (the ICCG set, the DuraPETase set, HotPETase positions) as an annotated overlay.
- Distance-to-triad and distance-to-cleft from the predicted structure, so mutations can be filtered by whether they are plausibly active-site, second-shell or surface.

---

## 6. Structure layer

For the top N candidates (start at 500, tune to disk budget):

1. Predict structure offline. ESMFold for throughput; Boltz-2 via **BoltzMaker** where a co-folded PET oligomer or MHET ligand is wanted.
2. Superpose onto IsPETase (PDB 6EQE and relatives) and record RMSD, TM-score.
3. Extract active-site geometry: triad distances and angles, oxyanion hole geometry, cleft width and depth, aromatic residue positions lining the cleft, solvent accessibility.
4. Store per-candidate mmCIF plus a geometry row in SQLite.

Cleft width and the aromatic clamp are the features most likely to separate real polyesterases from soluble-ester-only esterases, so surface them prominently.

---

## 7. Web app specification

House style throughout: marcdeller.com brand theme, Inter and Roboto Mono, `--md-primary` `#1e73be`, mandatory brand header, mobile responsive, British English, no em dashes.

Plotly for all charts (responsive config, transparent backgrounds). **Mol\*** for all 3D, loaded from CDN, one shared viewer component reused across tabs.

### Tabs

**Home**
What PANTS is, the therapeutic framing, headline counts (assemblies scanned, candidates recovered, candidates scored, candidates with predicted structures), and the current top 10 by therapeutic composite. One hero Plotly plot: activity probability versus predicted Topt, coloured by source environment, with the 37 °C therapeutic window shaded.

**Catalogue** (`/catalogue`)
The main table. Sortable, filterable, paginated, Tabulator-style. Columns: candidate ID, source environment, taxonomy, length, activity probability, predicted Topt, expressibility, aggregation, therapeutic composite, identity to nearest characterised enzyme, structure available.
Live-updating composite weight sliders in a sidebar. Filters for environment, taxon, identity band, triad completeness.

**Candidate** (`/candidate/<id>`)
Per-candidate detail. Mol\* viewer with the predicted structure, triad highlighted, cleft surface toggle, ligand overlay where a Boltz-2 co-fold exists. Alongside it: per-residue Plotly track (conservation, pLDDT, distance to triad), the score panel with each head's value and calibration, sequence with triad and oxyanion hole annotated, and nearest characterised relatives with identity.

**Compare** (`/compare`)
Two to four candidates side by side. Superposed Mol\* view, radar or parallel-coordinates Plotly of the score axes, aligned sequence view of the cleft-lining positions.

**Engineer** (`/engineer`)
Pick a candidate or paste a sequence. Plotly heatmap of variant effect scores (position by amino acid), with catalytic positions masked and literature-known beneficial substitutions marked. Mol\* view with the selected mutation highlighted in context. Variant basket: accumulate mutations, see the combined predicted effect, export a FASTA of the designed construct.

**Submit** (`/submit`)
Paste or upload a sequence. Runs the recall check and the ESM-2 35M heads live. Returns scores plus percentile against the catalogue, nearest catalogue neighbours, and nearest characterised enzymes. Under 10 seconds on CPU.

**Analysis** (`/analysis`)
Cross-catalogue views: embedding UMAP or PCA coloured by predicted activity with characterised enzymes overlaid as anchors, environment enrichment, taxonomy breakdown, the industrial-versus-therapeutic Topt distribution shift, and the retrieval baseline comparison (model rank versus E-value rank on held-out characterised enzymes).

**Methods** (`/methods`)
Full transparency page: data sources and versions, negative set construction, PU class prior, model architecture, calibration curves, evaluation protocol, and an explicit limitations section lifted from section 9 below. Non-negotiable given how easily this kind of tool is over-read.

---

## 8. Evaluation protocol

- **Splits by sequence cluster** (MMseqs2 at 30% and 50% identity), never by sequence. Report both.
- **Leave-one-family-out** across ESTHER families, to test whether the model generalises beyond the families it has seen.
- **Retrieval baseline**: model rank versus HMMER E-value rank on held-out characterised enzymes. This is the bar.
- **Calibration**: reliability diagrams, not just AUC. A calibrated probability is the point of the exercise.
- **Prospective set**: hold out any PETase characterised after a fixed date as a true blind test, in the same spirit as the AlphaFraud cutoff logic.
- Report performance separately for "quantitatively characterised" and "binary annotation only" subsets.

---

## 9. Limitations to state explicitly in the app

1. Positives number in the low hundreds. Every score is an extrapolation from a small, biased sample.
2. Published activity data is not harmonised across assay formats, and absolute rate predictions should not be trusted.
3. Crystalline PET degradation at 37 °C by any known enzyme is slow. PANTS ranks relative promise, not therapeutic viability.
4. Predicted structures are predictions. Cleft geometry from ESMFold on a metagenomic sequence with no close homologue carries real uncertainty.
5. Nothing here addresses delivery, immunogenicity, biodistribution, or what happens to liberated TPA and EG in vivo. Those decide whether any of this is a therapy.
6. Metagenomic candidates may be from unculturable organisms, may not express in a standard host, and may be fragments or misassemblies.

---

## 10. Open decisions for the plan to resolve

1. Which metagenome collections to scan first, given droplet storage limits. Plastisphere and landfill first, or breadth first?
2. Whether MHETase candidates share the pipeline or get a parallel one with their own seed and negatives.
3. Class prior estimate for the PU setup, and how sensitive rankings are to it.
4. Whether to attempt a quantitative activity regression head at all, or ship ranking only in v1.
5. Structure budget: how many candidates get ESMFold, how many get Boltz-2 co-folds, and the disk cost of the mmCIF assets.
6. Whether `/submit` and `/engineer` are public or gated, given compute cost per request.
7. Licensing, and whether the candidate catalogue is released as a dataset alongside the app.

---

## 11. Repository structure

```
pants/
  pipeline/
    recall/          # HMM profiles, MMseqs2 and HMMER drivers
    negatives/       # ESTHER harvest, matching, hard negative assembly
    embed/           # ESM-2 batch embedding
    train/           # head training, PU loss, calibration
    structure/       # ESMFold / BoltzMaker drivers, geometry extraction
    harmonise/       # literature activity extraction and normalisation
    db/              # SQLite schema, migrations, manifest writing
  app/
    __init__.py      # Flask app factory
    views/           # one module per tab
    static/
      css/           # house theme
      js/            # Plotly config, Mol* wrapper, shared components
      structures/    # precomputed mmCIF assets
    templates/
  models/            # trained head weights, calibration objects
  data/
    raw/
    interim/
    processed/
  eval/              # splits, metrics, baseline comparison
  deploy/            # systemd units, nginx config, gunicorn config
  tests/
```

---

## 12. Stack and style constraints

- **Back end**: Python, Flask, gunicorn, SQLite, nginx, certbot, systemd timer for the batch pipeline. Mirror AlphaFraud's deployment exactly.
- **Scientific Python**: biotite or gemmi for structure handling, HuggingFace transformers for ESM-2, scikit-learn for heads and calibration, numpy and scipy.
- **Front end**: server-rendered templates plus vanilla ES6. Plotly.js for charts, Mol\* for 3D, Tabulator for tables, Papa Parse for any CSV upload. **No React, Vue, npm, webpack, Streamlit or Dash.**
- **Design**: marcdeller.com theme, mandatory brand header with both links, mobile responsive, sticky header that does not overlap content, axis labels with units, threshold lines labelled.
- **Language**: British English throughout. No em dashes, use colons or parentheses.
- Every pipeline stage writes a manifest with input hashes, tool versions and model version. The app displays the current data version in the footer.

---

## 13. First task for the planning pass

Produce an implementation plan for a **v1 that is browsable end to end on a small dataset**: one metagenome collection, the recall stage, the hard negative set, one trained activity head with calibration, structures for the top 50 candidates, and the Home, Catalogue, Candidate and Methods tabs live at `pants.mdeller.com`. Engineer, Compare, Submit and Analysis come in v2.

Include a task breakdown with dependency order, the concrete Python packages per stage, the SQLite schema, an estimate of offline compute time and disk footprint, and an explicit check that on-demand inference fits inside the droplet's memory alongside gunicorn. Flag anything in section 10 that blocks v1 as opposed to blocking later versions.
