# 🧬 PANTS

> **Find PET-degrading enzymes that work at 37 °C in serum, not at 70 °C in a reactor.**

![python](https://img.shields.io/badge/python-3.14-3776AB?logo=python&logoColor=white) ![flask](https://img.shields.io/badge/flask-3.1-000000?logo=flask&logoColor=white) ![sqlite](https://img.shields.io/badge/sqlite-WAL-003B57?logo=sqlite&logoColor=white) ![esm2](https://img.shields.io/badge/ESM--2-t12--35M-467FF7) ![status](https://img.shields.io/badge/status-in%20development-fcb900) ![author](https://img.shields.io/badge/author-Marc%20C.%20Deller%2C%20D.Phil.-1C244B)

<table>
<tr>
<td>🌐 <b>Website</b></td><td><a href="https://marcdeller.com" target="_blank" rel="noopener noreferrer">marcdeller.com</a></td>
<td>✉️ <b>Contact</b></td><td><a href="mailto:marc@marcdeller.com">marc@marcdeller.com</a></td>
<td>🐙 <b>GitHub</b></td><td><a href="https://github.com/bellcheddar/PANTS" target="_blank" rel="noopener noreferrer">bellcheddar/PANTS</a></td>
</tr>
</table>

---

PANTS (PETase ANnotation and Triage System) discovers, triages and engineers PETase-like and MHETase-like enzymes for **therapeutic** use: degrading PET microplastics under physiological conditions rather than in an industrial reactor. It mines metagenomic sequence space, ranks candidates on axes that matter at 37 °C, and will serve the result as an interactive catalogue with predicted structures and in-silico engineering.

**Why it matters:** essentially the entire published PETase field optimises for industrial conditions (above PET's glass transition at roughly 70 °C, often alkaline, with thermostability and enzyme cost as the dominant concerns). Therapeutic use inverts almost all of it: peak activity at 37 °C, pH 7.4, in serum, against highly crystalline aged microplastic, with immunogenicity and protease resistance suddenly central. A naive homology search seeded on characterised PETases ranks candidates *towards* the industrial optimum, because that is where the well-annotated, heavily-published enzymes sit. Correcting for that bias is most of what PANTS is for. It is useful for: anyone triaging polyesterase sequence space on physiological rather than industrial criteria, and anyone who wants a calibrated probability rather than an E-value rank before committing a wet-lab screening budget.

## 🧭 Why this is not a homology search

Detection is solved. Discrimination is not.

A profile HMM built from characterised PET hydrolases returns thousands of α/β-hydrolase fold members sharing the Ser-His-Asp triad and the oxyanion hole. Almost none have meaningful activity on crystalline PET, and sequence identity to IsPETase is a weak predictor of PET activity, so E-value rank is close to uninformative about the property of interest.

The architecture that follows: **retrieval is the recall stage, the learned model is the precision stage.** MMseqs2 and HMMER cast the net exhaustively and interpretably; the model ranks within it on axes retrieval is blind to. That also yields a free baseline: if the model cannot beat E-value rank on held-out characterised enzymes, that is a cheap and publishable negative result.

## ⚗️ The therapeutic constraint

| Axis | Industrial optimum | Therapeutic requirement |
|---|---|---|
| Temperature | 65 to 75 °C | Peak activity at 37 °C |
| pH | 8 to 9 | 7.2 to 7.4 |
| Substrate | Amorphous, pre-treated, high surface area | Highly crystalline aged microplastic and nanoplastic |
| Medium | Buffer | Serum, plasma proteins, lipids, physiological ionic strength |
| Stability concern | Thermal denaturation over days | Serum protease resistance, aggregation, clearance |
| Immunogenicity | Irrelevant | Central |
| Product handling | Recovered and recycled | TPA and EG must be tolerable at achievable local concentrations |

## 🚧 Current status

**Pre-deployment.** The offline pipeline runs end to end from metagenome FASTA to embedded candidates. Nothing is live at `pants.mdeller.com` yet.

| Phase | State |
|---|---|
| 0: Scaffold, schema, manifest provenance | ✅ Complete |
| 1: Curation, hard negatives, activity data | ✅ Complete, gate still MARGINAL |
| 2: Recall (profile HMMs, MMseqs2, HMMER) | ✅ Complete |
| 4: ESM-2 embedding | ✅ Complete |
| 5: Activity head, calibration, evaluation | ⬜ Next |
| 6: Structures and active-site geometry | ⬜ Not started |
| 7 to 8: Web app and deployment | ⬜ Not started |

What is in the database today:

| Set | Count | Notes |
|---|---|---|
| **Candidates** | **128** | Mined from 2.2M metagenomic proteins, all triad-complete |
| Positives | 529 | Of which **16 experimentally evidenced**, the rest predicted (see below) |
| Hard negatives | 131 | Matched on five axes |
| Near misses | 125 | ESTHER `Cutinase` family: the decision boundary |
| Activity measurements | 47 | Km, Topt, pH optimum, each citing its PubMed IDs |
| Embeddings | 848 | ESM-2 t12-35M, 480-dim, frozen |
| Excluded from training | 70 | Fragments and length outliers, marked not deleted |

### Positives by evidence tier

The count that matters is not 529 but **16**: the number with experimental evidence behind the label.

| Tier | n | What it means |
|---|---|---|
| `EC-auto-annotated` | 449 | EC 3.1.1.101 assigned by similarity (ECO:0000256). A prediction, not a measurement |
| `ESTHER-family-predicted` | 50 | Family membership only |
| `ESTHER-family-protein-evidence` | 14 | Family, protein observed |
| `EC-experimental` | 10 | EC 3.1.1.101 with ECO:0000269 and PubMed citations |
| Curated wild types and variants | 6 | Hand-curated, sequence-verified, mutations validated |

## 🧪 What Phase 1 found

Two findings that changed the plan, both surfaced by the pre-training sanity gate rather than by review. Full detail in [`PHASE1_FINDINGS.md`](PHASE1_FINDINGS.md).

**The nine curated positives are one cluster, not nine examples.** They collapse into a single cluster at both 30% and 50% identity, because the engineered variants are 97.6% to 99.3% identical to their parents. Under the project's own evaluation rule (split by cluster, never by sequence) that is one independent example, so no cluster-split evaluation is possible over them and any cross-validation is pure leakage. The first ungrouped trivial baseline returned AUC 0.9996 for exactly that reason. Harvesting the full ESTHER polyesterase family took the positive set to 87 sequences in 11 clusters at 30%, at the cost of those being annotation-only labels.

**The hard negatives were separable on amino-acid composition alone.** Cluster-grouped, a classifier using nothing but 20 amino-acid fractions and length scored AUC 0.954 against a null of 0.495. The coefficients diagnosed it: negatives Leu-rich, positives Ser/Thr/Gly/Pro-rich, which is a **secreted-versus-cytoplasmic** signature rather than polyester chemistry. Every characterised polyesterase is secreted; the negative families were largely intracellular.

Negatives are now matched on five axes: length distribution, identity to nearest positive, genus cap, **signal peptide**, and **phylum**. Phylum matching specifically contributed almost nothing (0.845 to 0.842), which is itself informative: the residual is not GC-driven.

**Curating real activity data moved it further.** `EC 3.1.1.101` is poly(ethylene terephthalate) hydrolase, a curator's assignment of measured function rather than a family guess, and harvesting it took the positive set from 87 sequences in 11 clusters to 529 in 29:

| Positive set | Clusters | Composition baseline | Verdict |
|---|---|---|---|
| Curated only | 1 | 0.9996 (leakage, not a measurement) | invalid |
| Plus ESTHER family | 11 | 0.842 | MARGINAL |
| Plus EC 3.1.1.101 | 29 | **0.778** | MARGINAL |

Still short of the 0.75 pass mark, but the trend confirms the diagnosis: much of the apparent shortcut was a small-sample artefact that shrinks as real diversity arrives.

Consequently the composition baseline is a **permanently reported metric** alongside the retrieval baseline, not merely a pre-training gate. Any claim the model makes has to clear both.

**A caution about the evidence tiers.** Of the 449 entries carrying EC 3.1.1.101 by automatic annotation, none is a measurement: they hold `ECO:0000256` (by similarity), not `ECO:0000269` (experimental). They were briefly labelled "unreviewed", which reads as a curation backlog rather than the substantive difference it is. Sixteen positives have experimental evidence. That is the number the Methods tab will report.

## 🔭 The profile library and the recall run

Recall is a two-stage funnel: MMseqs2 casts the net across millions of sequences fast, HMMER makes the sensitive call on the survivors. Every candidate keeps its retrieval numbers (E-value, bitscore, profile identity), because those are the baseline the learned model has to beat.

The library is **one profile HMM per 30% sequence cluster**, each with its own catalytic anchor, rather than a single pooled profile. That matters: a single profile built over the polyesterases scored **0 of 111 near misses** as triad-complete, not because classic cutinases lack a catalytic triad but because they never aligned well enough for the columns to map. Per-cluster profiles took that to 79%, so the near misses survive recall and reach the scoring stage where they belong.

Anchors come from UniProt's own `Active site` annotation rather than being hardcoded. Cross-checked before adoption: aligning to a pooled profile and reading IsPETase's verified S160/D206/H237 columns predicted LCC as S165/D210/H242 and TfCut2 as S170/D216/H248, and UniProt's independently curated annotations give exactly those numbers.

| | |
|---|---|
| Library built from | 529 positives in 29 clusters at 30% identity |
| Profiles | 3 (264, 73 and 3 sequences), anchored on LCC, `P9WP41` and `A6WFI5` |
| Proteins scanned | 2,220,462 |
| Candidates recovered | **128**, all triad-complete |
| Runtime | 1,424 s (24 min) on an M1 Max |

### Funnel

| Stage | Surviving | Note |
|---|---|---|
| Scanned | 2,220,462 | |
| MMseqs2 prefilter, E ≤ 1e-5 | 1,698 | 0.08% |
| Matched a profile (hmmscan) | 892 | |
| Complete catalytic triad | 134 | 15% of profile-matched |
| Unique candidates written | **128** | content-addressed, so the same protein found in two assemblies collapses to one row |

Retention is 0.006%. That is a **choice** (a strict E-value against 500 positives, then a hard triad requirement), not a property of the data.

### By source environment

| Environment | Study | Proteins scanned | Candidates | Per 1M | Median length | Median identity | Max identity | Best bitscore |
|---|---|---|---|---|---|---|---|---|
| Compost | MGYS00006036, MGYS00005026 | 1,020,575 | 69 | 67.6 | 287 aa | 0.341 | 1.000 | 486.0 |
| Marine plastisphere | MGYS00006544 | 737,027 | 44 | 59.7 | 246 aa | 0.265 | 0.397 | 64.9 |
| Landfill | MGYS00004882 | 436,229 | 15 | 34.4 | 294 aa | 0.328 | 0.777 | 297.8 |
| Wastewater | MGYS00004904 | 26,631 | 0 | 0.0 | | | | |

**The identity bands are the interesting part**, because they separate rediscovery from genuinely unexplored sequence space:

| Environment | ≥70% identity (rediscovery) | 40 to 70% | <40% (novel) |
|---|---|---|---|
| Compost | 15 | 14 | 40 |
| Marine plastisphere | **0** | **0** | **44** |
| Landfill | 1 | 2 | 12 |

Compost gives the highest yield per million proteins and hands back 15 near-identical copies of enzymes that are already characterised, one of them a 100% identity match. That is unsurprising: LCC itself is leaf-branch compost derived, so compost is where the field has already looked.

**Every single marine plastisphere candidate sits below 40% identity to anything characterised**, with a best bitscore of 64.9 against compost's 486. Nothing in that cohort is a rediscovery. This is spec section 2's thesis in one table: E-value rank pushes the well-known enzymes to the top and the unexplored ones down, and re-ranking that is exactly what the learned model is for.

The wastewater assembly returned nothing, which is a reasonable null: it is the only source in the set that is neither plastic-associated nor compost.

## 🧱 Stack

| Layer | Choice |
|---|---|
| Recall | MMseqs2 18-8cc5c, HMMER 3.4 (brew binaries, shelled out to) |
| Embedding | ESM-2 `t12-35M`, frozen, CPU |
| Heads | scikit-learn logistic regression, PU-corrected, Platt/isotonic calibration |
| Structures | ESMFold offline, Boltz-2 via BoltzMaker for ligand co-folds (v2) |
| Storage | SQLite (WAL), one file, no external database |
| Web | Flask + gunicorn behind nginx, server-rendered templates plus vanilla ES6 |
| Front end | Plotly.js, Mol\*, Tabulator. No React, Vue, npm, webpack, Streamlit or Dash |

**Two virtual environments, deliberately.** The droplet has 3.8 GB of RAM shared with five other applications, so the always-on web process never imports torch: `requirements-web.txt` is Flask, gunicorn and gemmi, and nothing else. All heavy compute is precomputed offline on an M1 Max and shipped as SQLite rows plus static mmCIF.

## 🔧 Installation

```bash
git clone https://github.com/bellcheddar/PANTS.git
cd PANTS

# external tools
brew install hmmer mmseqs2

# pipeline venv (offline batch work: torch, transformers, scikit-learn)
python3 -m venv .venv
.venv/bin/pip install -r requirements-pipeline.txt

# web venv (always-on serving: no torch, ever)
python3 -m venv .venv-web
.venv-web/bin/pip install -r requirements-web.txt

# keep bulk data out of the iCloud-synced Documents tree
mkdir -p ~/PANTSData/{raw,interim}
ln -s ~/PANTSData/raw data/raw
ln -s ~/PANTSData/interim data/interim

cp .env.example .env
```

## 🚀 Usage

```bash
.venv/bin/python PANTS.py init                # create dirs, DB schema, check symlinks
.venv/bin/python PANTS.py curate-seeds        # fetch wild types, derive variants
.venv/bin/python PANTS.py harvest-negatives   # ESTHER hard negatives, matched
.venv/bin/python PANTS.py status              # database summary
.venv/bin/python PANTS.py serve               # local dev server on :8005
.venv/bin/pytest                              # 28 tests
```

| Command | Does |
|---|---|
| `init` | Creates the directory tree and the schema, enables WAL, warns if `data/raw` is not a symlink out of iCloud |
| `curate-seeds` | Fetches each wild type from UniProt by accession and derives every confirmed variant from its parent |
| `harvest-negatives` | Streams the ESTHER slice from UniProt, classifies by family, selects a matched negative set |
| `status` | Counts of candidates, scores, structures, characterised enzymes and recent runs |
| `serve` | Flask dev server (production uses gunicorn on port 8005) |

## 🔬 Sequences are fetched, never typed

Every sequence enters through the UniProt REST client. Engineered variants have no accession of their own, so each is stored as a **parent plus a mutation list** and the sequence is derived, with `apply_mutations` refusing any substitution whose stated parent residue does not match.

That check earns its keep: a wrong mutation set yields a sequence that is still a valid protein, still folds, still embeds and still trains. The error would never surface as a crash, only as quietly degraded scores. Where a complete mutation set could not be confirmed, the variant is recorded with **no sequence** and excluded from training, because a partial set gives a wrong sequence and that is worse than an honest gap.

`find_offset` determines the mature-versus-precursor numbering shift rather than guessing it: only an offset satisfying every mutation at once is accepted. All four confirmed variants validated at offset 0, across 14 residue positions.

Corrections this caught during curation: `Q6A0I4` was initially curated as Cut190 and is actually **TfCut2** (*Thermobifida fusca*). Cut190 is `W0TJ64`, and its strain assignment (AHK190 versus type strain P101, both 304 aa) is still unresolved.

## 📊 Evaluation protocol

| Element | Rule |
|---|---|
| Splits | By sequence cluster at 30% and 50% identity, never by sequence. Both reported |
| Generalisation | Leave-one-family-out across ESTHER families |
| Retrieval baseline | Model rank versus HMMER E-value rank on held-out characterised enzymes |
| Composition baseline | Amino-acid composition plus length, cluster-grouped. Reported permanently |
| Calibration | Reliability diagrams and Brier score, not just AUC |
| Prospective set | Any PETase characterised after a fixed date, held out as a blind test |
| Subsets | Reported separately for measured-activity and annotation-only positives |

## ⚠️ Limitations

1. Positives number in the low hundreds, and most carry family annotation rather than measured PET activity. Every score is an extrapolation from a small, biased sample.
2. Published activity data is not harmonised across assay formats. Absolute rate predictions should not be trusted.
3. Crystalline PET degradation at 37 °C by any known enzyme is slow. PANTS ranks relative promise, not therapeutic viability.
4. Predicted structures are predictions. Cleft geometry from ESMFold on a metagenomic sequence with no close homologue carries real uncertainty.
5. Nothing here addresses delivery, immunogenicity, biodistribution, or what happens to liberated TPA and EG in vivo. Those decide whether any of this is a therapy.
6. Metagenomic candidates may come from unculturable organisms, may not express in a standard host, and may be fragments or misassemblies.
7. The composition baseline sits at AUC 0.778. Until a model clears that as well as the E-value baseline, no claim of learned discrimination is supported.
8. Of 529 positives, only 16 carry experimental evidence. The rest are automatic EC annotation or family membership, so any head trained today is trained mostly on predicted labels.
9. Everything of interest is packed tightly in embedding space (characterised PET enzymes sit at cosine 0.96 or above to each other, candidates at a median 0.931 to their nearest known enzyme). The head discriminates small differences inside a dense cluster, not well-separated groups.

## 📚 Data sources

| Source | Use |
|---|---|
| UniProt / UniRef | Reference sequence space, taxonomy, evidence level, signal peptides |
| ESTHER | α/β-hydrolase family assignment, hard negatives, near misses |
| PAZy | Characterised plastic-degrading enzymes: the measured-activity positives |
| PDB | Experimental structures and ground-truth geometry |
| MGnify, JGI IMG/M, OceanDNA, Tara Oceans | Metagenomic assemblies for mining |
| Meltome Atlas, FireProtDB | Thermostability transfer learning |
| AlphaFold DB | Precomputed structures where a UniProt match exists |

## ✅ To Do

Roadmap for PANTS, roughly in dependency order. Suggestions welcome.

- [x] **Repository scaffold and the two-venv split.** `requirements-web.txt` carries Flask, gunicorn and gemmi and nothing else, so the always-on droplet process never imports torch. The droplet has 3.8 GB shared with five other applications, which makes this a memory constraint rather than a style preference
- [x] **SQLite schema with manifest provenance.** Thirteen tables, WAL enabled once at creation, every pipeline stage opening a run and writing input/output hashes, tool versions, git commit and wall time to both a table and a JSON file. A stage that raises still leaves its manifest behind
- [x] **Verify torch on Python 3.14 (plan risk 8).** torch 2.13.0 installs cleanly with MPS available, so no fallback to 3.11 was needed
- [x] **Install and pin the external tools.** HMMER 3.4 and MMseqs2 18-8cc5c, shelled out to rather than bound as libraries, with versions captured in every manifest
- [x] **Move bulk data out of the iCloud tree.** `data/raw` and `data/interim` symlink to `~/PANTSData`; macOS "Optimize Mac Storage" evicts large files mid-run and this machine has only ~62 GB free
- [x] **Curate the characterised seed set.** Wild types fetched from UniProt by accession; engineered variants derived from parent plus mutation list, with every substitution checked against the parent residue it names. All four confirmed variants validated at offset 0 across 14 positions
- [x] **Harvest ESTHER hard negatives.** Matched on five axes: length distribution, identity to nearest positive, genus cap, signal peptide and phylum
- [x] **Run the trivial-baseline gate before any embedding work (plan risk 1).** It fired, and found both that the curated positives are one cluster rather than nine examples, and that the negatives were separable on a secreted-versus-cytoplasmic composition signature
- [x] **Record an evidence level on every positive.** UniProt `protein_existence` separates the 23 with protein-level evidence from the 56 predicted or inferred, so the two are never pooled in a reported metric
- [x] **Make the composition baseline a permanent reported metric.** Stored alongside the retrieval baseline in `training_runs`, with `n_positive_clusters` recording independent units rather than the raw count
- [x] **Curate measured activity data.** Taken from UniProt's machine-readable, citable annotations rather than transcribed from PDFs, which is where fabrication risk lives. `EC 3.1.1.101` (poly(ethylene terephthalate) hydrolase) gave 459 entries and took the positive set from 87 sequences in 11 clusters to 529 in 29. 47 measurements extracted (21 Km, 8 Topt, 8 pH optima, 10 qualitative), each carrying its PubMed IDs, with `comparable_group_id` keyed on parameter plus substrate so a Km on pNP-butanoate is never pooled with one on PET film
- [ ] **Extend curation to rates on PET itself.** The extracted Km values are on soluble ester proxies (pNP esters), not on PET film or powder. Rates on real PET remain locked in paper supplementaries, and spec section 5.4's ordinal-within-paper fallback has not been built
- [ ] **Resolve the Cut190 strain ambiguity.** `W0TJ64` versus `C7MVE8`, both 304 aa, AHK190 versus type strain P101. Currently flagged in the seed notes and unresolved
- [ ] **Confirm the outstanding mutation sets.** DuraPETase, HotPETase, TurboPETase, Z1-PETase and Cut190\*\*SS are recorded without sequences because their complete mutation sets were not confirmed. A partial set yields a wrong sequence, which is worse than an honest gap
- [x] **Choose and acquire the metagenome collections, size-checked first.** 2,220,462 predicted proteins (858 MB) from landfill, marine plastisphere and compost assemblies. Only assemblies carry proteins: MGnify's largest plastisphere study has 357 samples and no protein sequences at all, being 16S amplicon
- [x] **Build the recall stage.** One profile HMM per 30% cluster, each anchored on UniProt's own Active site annotation, MMseqs2 prefilter then hmmscan and a triad completeness filter. 128 candidates from 2.2M proteins in 24 minutes, with discard counts reported at every step
- [ ] **Detect the oxyanion hole properly.** Currently a weak sequence proxy: the hole is formed by backbone amides, which is a structural property, so the real determination has to wait for the structure stage
- [x] **Embed the candidate set.** ESM-2 t12-35M, frozen, CPU, mean-pooled with padding and CLS/EOS excluded. 848 vectors at 480 dimensions in under a minute
- [x] **Filter fragments and length outliers.** From UniProt's own Fragment flag rather than a length cutoff, plus a 200 to 450 aa window derived from the experimentally evidenced positives. Marked, never deleted, so the catalogue stays complete and the exclusion stays auditable
- [ ] **Train the PET activity head.** PU-corrected loss with the class prior estimated rather than assumed, sensitivity-tested across 1/3/5/10%, then Platt or isotonic calibration
- [ ] **Run the full evaluation protocol.** Cluster splits at 30% and 50%, leave-one-family-out, retrieval baseline, reliability diagrams, prospective holdout by date, and separate reporting for the measured-activity and annotation-only subsets
- [ ] **Smoke-test ESMFold on one structure before committing to fifty.** Apple Silicon support is asserted rather than verified; Boltz-2 via BoltzMaker is the fallback, with real timings already measured on this machine
- [ ] **Extract active-site geometry.** Triad distances and angles, oxyanion hole, cleft width and depth, aromatic clamp and solvent accessibility, validated against known PDB structures before being run on any predicted one
- [ ] **Build the v1 tabs.** Home, Catalogue, Candidate and Methods, with Mol\* over static mmCIF and Plotly for every chart
- [ ] **Deploy to `pants.mdeller.com`.** Mirror AlphaFraud's gunicorn, nginx, systemd and certbot setup on port 8005, then add the entry to the mdeller.com launcher
- [ ] **Measure real droplet headroom before any v2 inference work.** The memory figures in the plan are estimates, not measurements
- [ ] **MHETase pipeline (v2).** Its own seed and negatives: MHETase is Tannase family, so a PETase-seeded profile search cannot reach it
- [ ] **Choose a licence**, and decide whether the candidate catalogue ships as a dataset alongside the application

## 📝 Licence

Not yet chosen. Licensing, and whether the candidate catalogue is released as a dataset alongside the application, is an open project decision.

---

## 👤 Author

**Marc C. Deller, D.Phil.**  
Structural biologist & drug discovery scientist  

<table>
<tr>
<td>🌐</td><td><a href="https://marcdeller.com" target="_blank" rel="noopener noreferrer">marcdeller.com</a></td>
<td>✉️</td><td><a href="mailto:marc@marcdeller.com">marc@marcdeller.com</a></td>
<td>🐙</td><td><a href="https://github.com/bellcheddar/PANTS" target="_blank" rel="noopener noreferrer">github.com/bellcheddar/PANTS</a></td>
</tr>
</table>
