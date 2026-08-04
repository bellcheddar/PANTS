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

**Pre-deployment.** Phase 0 (scaffold) is complete and Phase 1 (data acquisition) is partly complete. Nothing is live at `pants.mdeller.com` yet.

| Phase | State |
|---|---|
| 0: Scaffold, schema, manifest provenance | ✅ Complete |
| 1: Seed curation and hard negatives | 🟡 Partly complete, gate not cleared |
| 2: Recall (HMM profiles, MMseqs2, HMMER) | ⬜ Not started |
| 3 to 6: Embedding, training, structures | ⬜ Not started |
| 7 to 8: Web app and deployment | ⬜ Not started |

What is in the database today:

| Set | Count | Notes |
|---|---|---|
| Curated positives | 9 | 5 wild types fetched from UniProt, 4 variants derived and validated |
| Family positives | 79 | ESTHER `Polyesterase-lipase-cutinase`, annotation only |
| Hard negatives | 252 | Matched on five axes (see below) |
| Near misses | 111 | ESTHER `Cutinase` family: the decision boundary |
| Recorded without sequence | 5 | DuraPETase, HotPETase, TurboPETase, Z1-PETase, Cut190\*\*SS |

## 🧪 What Phase 1 found

Two findings that changed the plan, both surfaced by the pre-training sanity gate rather than by review. Full detail in [`PHASE1_FINDINGS.md`](PHASE1_FINDINGS.md).

**The nine curated positives are one cluster, not nine examples.** They collapse into a single cluster at both 30% and 50% identity, because the engineered variants are 97.6% to 99.3% identical to their parents. Under the project's own evaluation rule (split by cluster, never by sequence) that is one independent example, so no cluster-split evaluation is possible over them and any cross-validation is pure leakage. The first ungrouped trivial baseline returned AUC 0.9996 for exactly that reason. Harvesting the full ESTHER polyesterase family took the positive set to 87 sequences in 11 clusters at 30%, at the cost of those being annotation-only labels.

**The hard negatives were separable on amino-acid composition alone.** Cluster-grouped, a classifier using nothing but 20 amino-acid fractions and length scored AUC 0.954 against a null of 0.495. The coefficients diagnosed it: negatives Leu-rich, positives Ser/Thr/Gly/Pro-rich, which is a **secreted-versus-cytoplasmic** signature rather than polyester chemistry. Every characterised polyesterase is secreted; the negative families were largely intracellular.

Negatives are now matched on five axes: length distribution, identity to nearest positive, genus cap, **signal peptide**, and **phylum**. That moved the baseline from 0.954 to 0.842, which is still above the 0.75 pass threshold. Phylum matching specifically contributed almost nothing (0.845 to 0.842), which is itself informative: the residual is not GC-driven.

Consequently the composition baseline is now a **permanently reported metric** alongside the retrieval baseline, not merely a pre-training gate. Any claim the model makes has to clear both.

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
7. The composition baseline sits at AUC 0.842. Until a model clears that as well as the E-value baseline, no claim of learned discrimination is supported.

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
