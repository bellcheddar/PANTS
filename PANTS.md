# 🧬 PANTS

> **Find PET-degrading enzymes that work at 37 °C in serum, not at 70 °C in a reactor.**

![python](https://img.shields.io/badge/python-3.14-3776AB?logo=python&logoColor=white) ![flask](https://img.shields.io/badge/flask-3.1-000000?logo=flask&logoColor=white) ![sqlite](https://img.shields.io/badge/sqlite-WAL-003B57?logo=sqlite&logoColor=white) ![esm2](https://img.shields.io/badge/ESM--2-t12--35M-467FF7) ![status](https://img.shields.io/badge/status-in%20development-fcb900) ![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21807353-1C244B?logo=doi&logoColor=white) ![data](https://img.shields.io/badge/data-CC--BY--4.0-9b51e0) ![author](https://img.shields.io/badge/author-Marc%20C.%20Deller%2C%20D.Phil.-1C244B)

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

## 🔀 How it works

```
  SOURCE ENVIRONMENTS                      SEED SET
  ┌──────────────────────┐                 ┌───────────────────────────────────┐
  │ compost      1.0M    │                 │ 6 wild types  (UniProt + UniParc)  │
  │ marine       0.7M    │                 │ 9 variants    (parent + mutations) │
  │ landfill     0.4M    │                 │ 2 variants    (PDB constructs)     │
  │ wastewater   0.03M   │                 │ 5 HGMPs       (SciDB deposit)      │
  │ human gut   12.1M    │                 │ 449 EC 3.1.1.101 (annotation)      │
  └──────────┬───────────┘                 └─────────────────┬─────────────────┘
             │  predicted proteins                           │
             │                                               ▼
             │                              ┌──────────────────────────────┐
             │                              │ cluster at 30% identity      │
             │                              │ → 1 profile HMM per cluster  │
             │                              │ → anchor each on UniProt's   │
             │                              │   own ACT_SITE annotation    │
             │                              └───────────────┬──────────────┘
             ▼                                              │
  ╔═══════════════════════════════════════════════════╗     │
  ║  RECALL                                           ║ ◄───┘
  ║                                                   ║
  ║   MMseqs2 prefilter        E ≤ 1e-5               ║   fast, exhaustive
  ║        │                                          ║
  ║        ▼                                          ║
  ║   hmmscan vs profile library                      ║   assigns each survivor
  ║        │                                          ║   to the family it
  ║        ▼                                          ║   actually resembles
  ║   triad filter    Ser·His·Asp connected in SPACE  ║   ← geometry, not motif
  ║        │                                          ║
  ╚════════╪══════════════════════════════════════════╝
           │  retains ~0.006%  ·  keeps E-value + bitscore
           ▼                     as the baseline to beat
     ┌───────────────┐
     │  CANDIDATES   │  content-addressed on the sequence,
     │               │  so the same protein found twice is one row
     └───┬───────┬───┘
         │       │
         │       └──────────────────────────┐
         ▼                                  ▼
  ┌─────────────────────┐        ┌────────────────────────────────┐
  │ ESM-2 t12-35M       │        │ ESMFold  (≤ 450 aa)            │
  │ frozen, CPU         │        │   │                            │
  │ 480-dim, mean-pool  │        │   ▼                            │
  └──────────┬──────────┘        │ superpose onto IsPETase 6EQE   │
             │                   │   │  ← at WRITE time, so the   │
             ▼                   │   ▼    browser does none       │
  ┌─────────────────────┐        │ active-site geometry           │
  │ activity head       │        │   triad distances, cleft width,│
  │ PU-corrected,       │        │   aromatic clamp               │
  │ calibrated          │        └────────────────┬───────────────┘
  └──────────┬──────────┘                         │
             │        ┌───────────────────────────┘
             ▼        ▼
        ┌────────────────────────────────────────────┐
        │  SQLite  →  Flask  →  Catalogue · Superpose │
        │                       Candidate · Methods   │
        └────────────────────────────────────────────┘

  Every stage writes a manifest: input hashes, tool versions, git commit, wall time.
```

Two things in that diagram are load-bearing and easy to miss. The **triad filter tests
whether Ser, His and Asp are connected in space**, not whether a motif matches, so three
residues present in sequence but not in contact correctly fail. And **superposition
happens when a structure is written**, not in the browser, which is what makes an
N-structure overlay cost N file loads instead of N alignments.

## 🌡️ The therapeutic gap, measured

The table above is the premise. These are the numbers. Every enzyme named in the project
brief, plus every other PET hydrolase with a published optimum that could be found.

Rows are split by **how the number was obtained**, because that difference matters more
than the value: the first block is extracted programmatically from UniProt's curated
`BIOPHYSICOCHEMICAL PROPERTIES` with `ECO:0000269` experimental evidence and PubMed IDs
attached, and is what sits in the database. The second is from the literature and is
recorded here for context only.

### Measured, in the database, each with its citation

| Enzyme | Topt | pH opt | Source |
|---|---|---|---|
| **HGMP01** (human gut) | **40 °C** | **~7.4, broad across pH 7.x** | PMID 39551294 |
| IsPETase | 40 °C | 9.0 | PMID 26965627, 29603535 |
| LCC | 50 °C | 8.5 | PMID 22194294, 24593046 |
| *T. alba* est1 (`D4Q9N1`) | 50 °C | 6.0 | PMID 25910960 |
| *T. alba* est2 (`F7IX06`) | 50 °C | 6.0 | PMID 20393707 |
| TfCut1 (`G8GER6`) | 55 °C | 8.0 | PMID 23604968 |
| TfCut2 | 55 °C | 8.0 | PMID 15638529, 20816933 |
| *T. fusca* (`Q47RJ6`) | 60 °C | 8.0 | PMID 18658138, 20729325 |
| *T. fusca* (`Q47RJ7`) | 60 °C | 8.0 | PMID 18658138 |

### The engineered lineage, from the literature

Every variant named in the project brief. They are in the database at **weaker
provenance** than the block above: `extraction_confidence='review'` and `ECO:0000305`
(inferred by curator) rather than `ECO:0000269` (experimental, from a publication), because
the values were collated from a secondary review rather than read from each primary paper.

**Two of the five missing sequences were recovered from crystal structures**, which is a
stronger route than deriving them: a PDB SEQRES is the construct that was actually
expressed, crystallised and assayed, so nothing is applied or assumed. Each was verified by
aligning against its parent and checking the substitution count against the published one,
rather than trusting the name on the entry.

| Variant | Source | Length | Substitutions vs parent | Published |
|---|---|---|---|---|
| HotPETase | PDB `7QVH` | 272 aa | 21 vs IsPETase | ~21 |
| Cut190\*\*SS | PDB `7CEF` | 262 aa | 4 vs Cut190 (incl. S226P/R228S) | S226P/R228S |

Both are stored as the **mature construct**, so they are shorter than their precursor parents.

**Three still carry no sequence: DuraPETase, TurboPETase and Z1-PETase.** None has a
name-matched PDB entry, none is in PAZy, and their mutation sets are in primary
supplementaries. The PDB does hold generic "IsPETase variant V20/V22" structures that
*might* be among them, but assigning one on a guess is exactly the error that would
silently poison the profiles, so they stay empty.

Values collated in
[Engineering Evolution: How FAST-PETase and Other Variants Are Transforming Plastic Biodegradation](https://marcdeller.com/engineering-evolution-how-fast-petase-and-other-variants-are-transforming-plastic-biodegradation/).

| Variant | Topt | Parent | Notes |
|---|---|---|---|
| Z1-PETase | 30 °C | IsPETase | 13 mutations, two engineered disulfides; 40x expression yield |
| IsPETase (wild type) | 30 to 35 °C | native | Weak on crystalline PET |
| **DuraPETase** | **37 °C** | IsPETase | 10 mutations; +31 °C thermostability, ~300x activity |
| FAST-PETase | 50 °C | ThermoPETase | 38x activity; 33.8 mM monomers in 96 h |
| HotPETase | 60 to 65 °C | IsPETase | 21 mutations; melting temperature 82.5 °C |
| Cut190\*\*SS | 65 °C | actinomycete cutinase | Calcium-dependent conformational switching |
| TurboPETase | 65 to 68 °C | BhrPETase M6 | 98.2% depolymerisation at 200 g/kg in 8 h |
| LCC-ICCG | 65 to 72 °C | LCC | 1.3 g PET waste in 3 days from 1.25 mg enzyme |

**A discrepancy worth stating rather than smoothing over:** UniProt's curated value for
IsPETase is 40 °C, while the literature review gives 30 to 35 °C. Both are defensible and
they were measured on different substrates under different assays, which is exactly the
harmonisation problem the project brief raises. The database keeps the UniProt value
because it carries an experimental evidence code and a citation; this table shows both
rather than picking a winner.

**Read the two tables together and the shape of the field is obvious.** The entire
engineered lineage runs *away* from body temperature: FAST-PETase to 50 °C, HotPETase to
60 to 65, TurboPETase and LCC-ICCG to 65 to 72. That is rational, because they were
optimised for a reactor above PET's glass transition.

Only three enzymes sit anywhere near 37 °C, and two of them (Z1-PETase, DuraPETase) are
engineered variants of IsPETase that reach it by trading away the activity the others
gained. **HGMP01 is the only one that arrived there naturally**, and the only measured
enzyme in this project with a near-neutral pH optimum as well: everything else wants pH 8
to 9, and even IsPETase asks for 9.

That is the therapeutic gap, and it is not an argument: it is what the numbers say. It is
also the reason the human gut was added as a fourth source environment, since an enzyme
resident at 37 °C and pH 7.4 has already been selected under something close to the
target conditions.

### Why HGMP01 matters, and the trap in retrieving it

HGMP01 is a PET hydrolase [identified from the human gut microbiome](https://pubmed.ncbi.nlm.nih.gov/39551294/)
(*Int J Biol Macromol* 283, 2024). It hydrolyses PET nanoparticles, outperformed the four
other candidates in that study, and shares only about 5% identity with IsPETase. Its
homologues are distributed across 41 families and 94 genera of gut microbes.

It is in no public sequence database: zero hits across UniProt, NCBI protein and NCBI
nuccore. The sequences live in the authors' SciDB deposit.

**The paper and the deposit number the enzymes differently, and matching on the name
would have mislabelled all five.** The paper's HGMP01 to HGMP05 are the deposit's HGMP03,
04, 06, 07 and 08. The mapping is pinned by two independent facts:

| Paper | Length | Deposit | Length | Identifier | Homologues |
|---|---|---|---|---|---|
| **HGMP01** | 275 | HGMP03 | 275 | **`GUT_GENOME238302_00589`** | **697** |
| HGMP02 | 341 | HGMP04 | 341 | `GUT_GENOME243637_00613` | 131 |
| HGMP03 | 323 | HGMP06 | 323 | `GUT_GENOME137663_00143` | 96 |
| HGMP04 | 282 | HGMP07 | 282 | `GUT_GENOME171691_00743` | 1000 |
| HGMP05 | 321 | HGMP08 | 320 | `GUT_GENOME244370_00064` | 87 |

Every length is distinct, so the assignment is 1:1 with no ambiguity. Independently, the
paper states that the homologue search "identified a total of 697 putative HGMP01-like
enzymes", and the deposit's HGMP03 returned exactly 697 hits. Two unrelated numbers
agreeing is what makes this a determination rather than a guess.

Every sequence is length-checked against the paper's own table before being stored, and
that check immediately caught a transposition between the similar identifiers
`GUT_GENOME244370_00064` and `GUT_GENOME244064_00699`, which would otherwise have filed
the wrong protein as HGMP05 with nothing downstream to reveal it.

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
| **Candidates** | **439** | Mined from 14.8M metagenomic proteins across four environments, all triad-complete |
| Positives | 854 | Of which **341 experimentally measured**, the rest predicted (see below) |
| Hard negatives | 131 | Matched on five axes |
| Measured-set head | AUC **0.976** | 45 independent clusters, against a 0.829 composition baseline |
| Near misses | 150 | 124 ESTHER `Cutinase` family, plus **26 within-family negatives**: PAZy enzymes measured on another plastic that share a 30% cluster with a PET-active one |
| Activity measurements | 48 | Km, Topt, pH optimum, each citing its PubMed IDs |
| Embeddings | 848 | ESM-2 t12-35M, 480-dim, frozen |
| Excluded from training | 70 | Fragments and length outliers, marked not deleted |

### Positives by evidence tier

The count that matters is not 854 but **341**: the number with a measurement behind the label.

| Tier | n | What it means |
|---|---|---|
| `EC-auto-annotated` | 449 | EC 3.1.1.101 assigned by similarity (ECO:0000256). A prediction, not a measurement |
| `ESTHER-family-predicted` | 50 | Family membership only |
| `ESTHER-family-protein-evidence` | 14 | Family, protein observed |
| `EC-experimental` | 10 | EC 3.1.1.101 with ECO:0000269 and PubMed citations |
| Curated wild types and variants | 6 | Hand-curated, sequence-verified, mutations validated |
| **`PAZy-measured`** | **312** | In PAZy because activity was **measured** on a plastic and published, each with a DOI |
| `HGMP-measured` | 5 | Human gut PET hydrolases with measured activity (PMID 39551294) |

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

### The labels were the blocker, and PAZy resolved it

For most of this project a head trained on the catalogue's own labels scored **AUC 1.000**,
which was a failure rather than a success: 449 positives carried `EC 3.1.1.101` from
`ECO:0000256`, meaning the label was assigned **by sequence similarity**. A sequence model
reproducing that is close to tautological. Meanwhile the positives with real experimental
evidence numbered 17 and spanned 5 clusters: too few for the cluster-split protocol to run
at all, so the honest head could not be scored rather than scoring badly.

[PAZy](https://api.pazy.eu/api) inverts the inclusion criterion. An enzyme is listed
because activity was **measured** on a plastic and published, with the DOI attached.

| | Measured positives | Clusters at 30% | Clusters at 50% |
|---|---|---|---|
| Before | 17 | 5 | 7 |
| With PAZy | 333 | 51 | 73 |
| Plus the confirmed engineered variants | **341** | **50** | **72** |

Trained on the measured set alone (300 after length filtering, 45 clusters, 220 negatives
and near misses):

| Metric | Value |
|---|---|
| AUC, cluster-grouped | **0.976 ± 0.021** over 4 valid folds |
| Average precision | 0.987 |
| Brier score | 0.052 |
| Composition-only baseline | 0.829 |

That is a different kind of number from the 1.000. The labels are experimental rather than
similarity-derived, the split is across 45 independent units rather than impossible, and
the head clears the composition baseline by a real margin instead of tying with it.

**What it still does not show, tested rather than assumed.** Running the head against the
125 near misses gave **AUC 1.000**, which is a warning rather than a triumph. Every near
miss belongs to a single ESTHER family (`Cutinase`), and one homogeneous family separates
as a block: composition alone scores 0.910 on the same contrast. That is a **family-level**
question wearing the clothes of a functional one.

| Contrast | ESM-2 head | Composition only |
|---|---|---|
| vs distant α/β-hydrolase families | 0.979 | 0.887 |
| vs near misses | **1.000** | 0.910 |
| vs both | 0.976 | 0.829 |

The honest test needs PET-**inactive** members of the polyesterase family, measured and
published. Those barely exist, and the reason is structural: **databases record what
works.** PAZy lists enzymes because activity was found; nobody systematically publishes
"we assayed this polyesterase on PET and it did nothing". The negative class for the
question that matters is largely unwritten, which is publication bias rather than a
curation gap, and no further database mining will fix it.

### 🎯 The within-family test, run

There is one partial way round that, and it has now been run rather than argued about.
PAZy curates enzymes measured on **other** plastics too (PA, PUR, PLA, PBAT, PHA). Most are
different folds doing different jobs, and including a nylon amidase would be another easy
negative. But 26 of them **share a 30% cluster with a PET-active enzyme**, which puts them
inside the polyesterase family boundary. That is the closest thing to the missing negative
class that exists.

Running the same head against negatives of increasing closeness, everything else fixed:

| Negative regime | Positives | Negatives | Head AUC | Composition only |
|---|---|---|---|---|
| Near-miss (one ESTHER `Cutinase` family) | 305 | 110 | 0.997 ± 0.006 | 0.890 |
| Out-of-family (distant α/β-hydrolases) | 305 | 109 | 0.975 ± 0.022 | 0.906 |
| **Within-family** (measured on another plastic) | 305 | 26 | **0.850 ± 0.055** | 0.651 |
| **Within-family, mixed clusters only** | 152 | 26 | **0.493 ± 0.297** | 0.398 |

**Discrimination decays monotonically as the negatives get closer, and collapses to chance
at the boundary that matters.** The head beats composition in every regime except the last.
In the strictest one, where every negative has a PET-active enzyme inside its own cluster,
so family membership carries no information at all, it scores 0.493: a coin flip.

**That last row is underpowered, and saying so is part of the result.** One cluster holds
85% of the positives, so the cluster-grouped folds come out badly lopsided, and the
per-fold AUCs are 0.193, 0.714, 0.208 and 0.857. A mean of 0.493 across that spread cannot
distinguish "no signal" from "not enough data to detect one". What it does rule out is a
large, easily-learned signal, which the 0.976 headline might otherwise be read as implying.

Two caveats travel with these numbers permanently:

- **The negatives are weak.** PAZy records only positive substrate associations, so "PET
  not listed" conflates *inactive* with *never assayed*. Some fraction are false negatives,
  which pushes measured discrimination **down**: the within-family figures are a lower
  bound, not a point estimate.
- **The composition baseline is doing more work than it appears.** It reaches 0.906 against
  out-of-family negatives. Most of the headline 0.976 is therefore available from amino
  acid composition alone, and the head's real contribution is the margin above that, not
  the absolute number.

The practical conclusion is not that the approach fails, but that **the evaluation is
label-limited rather than method-limited**, and the binding constraint is now a few dozen
measured within-family negatives rather than any modelling choice.

Three routes that would, in [`PHASE1_FINDINGS.md`](PHASE1_FINDINGS.md): ordinal
within-paper rankings, regression on measured rates instead of a binary label, and
geometry as an annotation-independent axis.

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

| Stage | Surviving |
|---|---|
| Scanned | 14,778,289 |
| Passed the prefilter, profile scan and triad filter | **439** |

Retention is 0.003%. That is a **choice** (a strict E-value, then a hard requirement that
Ser, His and Asp be connected in space), not a property of the data.

### By source environment

| Environment | Proteins scanned | Candidates | Per 1M | Median %ID to nearest characterised |
|---|---|---|---|---|
| **Human gut** | **12,584,458** | **311** | 25 | 29.7% |
| Compost | 1,020,575 | 69 | 68 | 34.1% |
| Marine plastisphere | 737,027 | 44 | 60 | 26.5% |
| Landfill | 436,229 | 15 | 34 | 32.8% |
| Wastewater | 26,631 | 0 | 0 | |
| **Total** | **14,778,289** | **439** | 30 | |

The gut is now 71% of the catalogue by count and 85% by sequences scanned. Its *yield* is
the lowest of the four productive environments, which is the honest reading: polyesterases
are rarer in the gut than in compost, and the gut cohort is large because it was sampled
deeply, not because it is rich.

**The identity bands are the interesting part**, because they separate rediscovery from genuinely unexplored sequence space:

| Environment | ≥70% identity (rediscovery) | 40 to 70% | <40% (novel) |
|---|---|---|---|
| Compost | 15 | 14 | 40 |
| Marine plastisphere | **0** | **0** | **44** |
| Landfill | 1 | 2 | 12 |

Compost gives the highest yield per million proteins and hands back 15 near-identical copies of enzymes that are already characterised, one of them a 100% identity match. That is unsurprising: LCC itself is leaf-branch compost derived, so compost is where the field has already looked.

**Every single marine plastisphere candidate sits below 40% identity to anything characterised**, with a best bitscore of 64.9 against compost's 486. Nothing in that cohort is a rediscovery. This is spec section 2's thesis in one table: E-value rank pushes the well-known enzymes to the top and the unexplored ones down, and re-ranking that is exactly what the learned model is for.

The wastewater assembly returned nothing, which is a reasonable null: it is the only source in the set that is neither plastic-associated nor compost.

### The blind test: did recall find HGMP01-like enzymes unaided?

HGMP01 was **not** used to seed recall. Its sequence was unavailable when the gut scan was
designed, and once obtained it was deliberately left out of the profiles, because the paper
that describes it reports HGMP01-like genes as widely distributed across the gut
microbiome. Recovering them without the seed is a real test; seeding first would have made
it circular.

Against the five measured HGMPs:

| | |
|---|---|
| Gut candidates | 311 |
| With a detectable hit to any HGMP | **254 of 311 (82%)** |
| Whose nearest HGMP is HGMP01 | **33** |
| Best identity to HGMP01 | 38.3% |
| Median identity to nearest HGMP | 27.5% |

So the pipeline independently recovered a large cohort of gut proteins resembling
enzymes it had never been shown, and 33 of them sit closest to the one with measured
activity at 40 °C and neutral pH.

**What this does not show.** None of them *is* HGMP01: the best match is 38.3%, so these
are relatives rather than rediscoveries. That is consistent with the source paper, which
found 697 HGMP01-like enzymes at 20% identity or better across the whole UHGP-100
database, where this scan covered 50 assemblies. It is also the expected outcome: HGMP01
comes from a UHGG genome that these particular assemblies need not contain.

The useful claim is narrower than "PANTS found HGMP01" and more interesting than nothing:
**recall reaches the right neighbourhood of sequence space unaided**, in an environment
nobody had asked it about until the therapeutic framing pointed there.

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
7. The composition baseline sits at AUC 0.829 on the measured set. The head clears it (0.976), but a model must keep clearing it, and the retrieval baseline, to claim learned discrimination.
8. AUC 0.976 is against hard negatives from other α/β-hydrolase families. Ranking PET activity **within** the polyesterase family has now been tested against the 26 PAZy enzymes measured on other plastics that share a 30% cluster with a PET-active one: discrimination decays to **0.850** and, restricted to mixed clusters where family membership carries no information, to **0.493** (chance). That strictest test is underpowered, and the negatives are weak because "PET not listed" conflates inactive with never assayed, so it is a lower bound rather than a verdict. The evaluation is label-limited, not method-limited.
9. Old note, kept: the composition baseline was 0.778 on the earlier annotation-heavy positive set. Until a model clears that as well as the E-value baseline, no claim of learned discrimination is supported.
10. Of 854 positives, 341 carry a measurement. The rest are automatic EC annotation or family membership, so any head trained today is trained mostly on predicted labels.
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
- [x] **Add the gut microbiome as a fourth source environment.** The three current environments (compost, marine plastisphere, landfill) are all external. A human gut metagenome is the one that matters most for the therapeutic framing: an enzyme already resident at 37 °C, pH 7.4 and in a proteolytic environment has been selected under something close to the target conditions, rather than being asked to work far from its optimum. **HGMP01** is the concrete starting point: it is named in spec section 1, it is metagenome-derived rather than a variant of a characterised parent, and it is currently the one named enzyme the recall stage is expected to recover from sequence space instead of being seeded with. Done: **12,584,458 gut proteins from five MGnify studies**, registered as `human_gut` and now 71% of the catalogue. **HGMP01 itself could not be seeded**: it returns zero hits in UniProt and in NCBI protein/nuccore, and its paper (PMID 39551294) has no linked sequence records, so the sequence sits in supplementary material. That is not a blocker but an advantage, because the same paper reports HGMP01-like genes as widely distributed across the gut microbiome, so recall finding them **unaided** is a blind test rather than a circular one
- [x] **Obtained all five HGMP sequences.** No public database holds them (zero hits across UniProt, NCBI protein and nuccore), but the authors deposited them on SciDB. **The paper and the deposit use different numbering**: the paper's HGMP01 to HGMP05 are the deposit's HGMP03, 04, 06, 07 and 08, so matching on name would have mislabelled all five. The mapping is pinned by two independent facts, sequence length (all five distinct, so 1:1) and homologue count (the paper's "697 putative HGMP01-like enzymes" matches the deposit's HGMP03 exactly). **HGMP01 = `GUT_GENOME238302_00589`, 275 aa, optimum 40 °C at near-neutral pH**: the only measured PET hydrolase in this project whose optimum is close to physiological
- [x] **Built the ordinal-within-paper fallback, and established why it is the route rather than a fallback.** Before this the database held exactly **one** quantitative value measured on PET; all 21 Km values were on soluble pNP-ester proxies. Attempting the real thing showed why: **PET rates live in figure panels**, not in extractable text or structured deposits, so they cannot be harvested at scale. The brief anticipated this and proposed within-paper ordinal ranking, which now exists as a `parameter_type` with `ordinal_rank_in_paper`, seeded from PMID 39551294 (HGMP01 first of five on PET nanoparticles, the other four recorded as equal-second because the paper does not order them). Assay conditions are attached so the comparison is interpretable
- [ ] **Get more measured within-family negatives.** This is now the binding constraint on the whole evaluation, ahead of any modelling work. The within-family test runs on 26 PAZy enzymes measured on other plastics, of which only 7 clusters overlap the PET-active set, and one cluster holds 85% of the positives, so the strictest contrast cannot separate "no signal" from "too little data". A few dozen polyesterases assayed on PET with a reported negative result would settle it. Candidate routes: the ordinal-ranking pass below (papers that assay several enzymes under one protocol often report the inactive ones), and PAZy entries whose primary papers include unreported negative controls
- [ ] **Populate ordinal rankings across the PAZy literature.** 312 PAZy entries carry DOIs; papers that assay several enzymes under one protocol rank them, and that ranking is extractable by reading where an absolute rate is not. This is the only route to a training signal that does not need a negative class, which the near-miss finding showed is otherwise unavailable
- [x] **Resolved the Cut190 strain ambiguity.** `W0TJ64`, not `C7MVE8`. Length could never separate them (both 304 aa), but the crystal structures can: **4WFI, 4WFJ, 4WFK, 5ZNO, 5ZRQ and 5ZRR all cross-reference W0TJ64**, and C7MVE8 has no PDB entry at all. Structural evidence beats a name match, and the seed was right by luck rather than by evidence until now
- [x] **Confirm the outstanding mutation sets.** All five are resolved, and two further variants (DepoPETase, LCC-A2) were found in the same pass. HotPETase and Cut190\*\*SS came from crystal structures; DuraPETase, TurboPETase and Z1-PETase from verified mutation sets, all matching at offset 0. Every engineered variant now carries a sequence, and the whole seed set rebuilds from one command
- [x] **Choose and acquire the metagenome collections, size-checked first.** 2,220,462 predicted proteins (858 MB) from landfill, marine plastisphere and compost assemblies. Only assemblies carry proteins: MGnify's largest plastisphere study has 357 samples and no protein sequences at all, being 16S amplicon
- [x] **Build the recall stage.** One profile HMM per 30% cluster, each anchored on UniProt's own Active site annotation, MMseqs2 prefilter then hmmscan and a triad completeness filter. 128 candidates from 2.2M proteins in 24 minutes, with discard counts reported at every step
- [ ] **Detect the oxyanion hole properly.** Currently a weak sequence proxy: the hole is formed by backbone amides, which is a structural property, so the real determination has to wait for the structure stage
- [x] **Embed the candidate set.** ESM-2 t12-35M, frozen, CPU, mean-pooled with padding and CLS/EOS excluded. 848 vectors at 480 dimensions in under a minute
- [x] **Filter fragments and length outliers.** From UniProt's own Fragment flag rather than a length cutoff, plus a 200 to 450 aa window derived from the experimentally evidenced positives. Marked, never deleted, so the catalogue stays complete and the exclusion stays auditable
- [ ] **Train the PET activity head.** PU-corrected loss with the class prior estimated rather than assumed, sensitivity-tested across 1/3/5/10%, then Platt or isotonic calibration
- [ ] **Run the full evaluation protocol.** Cluster splits at 30% and 50%, leave-one-family-out, retrieval baseline, reliability diagrams, prospective holdout by date, and separate reporting for the measured-activity and annotation-only subsets
- [x] **Smoke-test ESMFold before committing to a full run.** Cleared: 8.44 GB model (it bundles ESM-2 3B), 31 min to load, 109 s to fold 290 residues. No Boltz-2 fallback needed. The prediction reproduces IsPETase's crystal active site exactly (S160/D206/H237, Ser-His 2.98 A against 2.94)
- [x] **Extract active-site geometry, validated on crystal structures first.** The triad is found by geometry rather than sequence position, and recovers the published triads of 6EQE and 4EB0 exactly. Validation caught a brittle cleft-width metric: a hard radius feeding a max meant a 0.4 A shift in one residue halved the answer, scoring IsPETase's own prediction at 11.09 A against its crystal's 20.90 A. Now a percentile over a wider radius, and crystal-to-prediction disagreement fell from 9.81 A to 1.81 A
- [ ] **Establish whether geometry actually tracks PET activity.** The robust metric separates the fungal cutinase from the Polyesterase-lipase-cutinase family by ~7 A, but the family members do not order among themselves by PET activity: IsPETase, the best ambient-temperature degrader of the three, has the narrowest cleft of them. Family-level separation is a much easier task than the one that matters, and only the former is currently demonstrated
- [x] **Build the v1 tabs.** Home, Catalogue, Candidate, Compare and Methods, with 3Dmol.js over static PDB and Plotly for every chart
- [x] **Superposed interactive structure viewer, first working version.** Verified in a real browser, not assumed: multiple structures load into one Mol\* viewport in distinct colours and genuinely overlay, which confirms the pre-superposed-at-write-time design end to end. Two rendering defects remain, both recorded below
- [x] **Render structures as cartoon rather than wireframe.** Root cause found: Mol\* fell back to lines because gemmi's mmCIF carries no polymer or secondary-structure annotation, and the UMD build exports only 9 symbols. Moved to 3Dmol.js and now write HELIX/SHEET records computed with biotite's P-SEA `annotate_sse`, which 3Dmol parses to set `atom.ss` directly
- [x] **Draw the catalytic triad.** Resolved by the same move: `molstar.MolScriptBuilder` is not exported by the UMD viewer build, so the selection expression could never be constructed. 3Dmol selects the triad residues directly
- [ ] **Design the viewer interaction properly.** Current version is deliberately minimal: pick structures, toggle the reference, reset the camera. How the overlay set is chosen, how the triad is surfaced, whether IsPETase shows as a ghost and whether pLDDT colouring becomes the default are all still open View candidate structures overlaid, interactively, with the catalytic triad highlighted. Promoted to a first-class deliverable because geometry is the only annotation-independent signal here: a user judging a candidate needs to look at it against IsPETase rather than trust a score. Structures are already superposed onto 6EQE at write time so the browser does no alignment, and triad residue numbers are stored per candidate. The interaction itself (how the overlay set is chosen, how the triad is surfaced, whether IsPETase shows as a ghost, whether pLDDT colouring is the default) is deliberately still to be designed
- [ ] **Deploy to `pants.mdeller.com`.** Mirror AlphaFraud's gunicorn, nginx, systemd and certbot setup on port 8005, then add the entry to the mdeller.com launcher
- [ ] **Measure real droplet headroom before any v2 inference work.** The memory figures in the plan are estimates, not measurements
- [ ] **MHETase pipeline (v2).** Its own seed and negatives: MHETase is Tannase family, so a PETase-seeded profile search cannot reach it
- [x] **Published the catalogue as a citable dataset.** Zenodo, CC BY 4.0, concept DOI [10.5281/zenodo.21807353](https://doi.org/10.5281/zenodo.21807353), nine files each size-verified on upload. Code licence still to choose

## 📦 Dataset

The catalogue is deposited on Zenodo under CC BY 4.0: candidates with their retrieval
scores and active-site geometry, the reference set split by evidence tier, the measured
activity data with its citations, every pipeline run with its discard counts, and the
predicted structures.

[![DOI](https://img.shields.io/badge/DOI-10.5281%2Fzenodo.21807353-1C244B?logo=doi&logoColor=white)](https://doi.org/10.5281/zenodo.21807353)

| | |
|---|---|
| Concept DOI (always latest) | [10.5281/zenodo.21807353](https://doi.org/10.5281/zenodo.21807353) |
| This version (v0.1.0) | [10.5281/zenodo.21807354](https://doi.org/10.5281/zenodo.21807354) |

Every figure in the deposit's description and datasheet was **recomputed from the
deposited files** by `scripts/build_release.py` rather than copied from this README.
That is deliberate: on a sibling project six figures reached a Zenodo description by
being copied from prose describing a checkpoint that was measured but never shipped, and
a DOI would have made them permanent.

**Cite as:** Deller, M. C. (2026). *PANTS: a triaged catalogue of candidate PET-degrading
enzymes for therapeutic conditions*. Version 0.1.0. Zenodo. https://doi.org/10.5281/zenodo.21807353

## 📝 Licence

**Data:** CC BY 4.0, via the Zenodo deposit above.

**Code:** not yet chosen. Attributing PANTS does not discharge the obligation to the
sources it derives from: UniProt and UniRef (CC BY 4.0), the PDB (CC0), PAZy, ESTHER and
MGnify, each of which carries its own terms.

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
