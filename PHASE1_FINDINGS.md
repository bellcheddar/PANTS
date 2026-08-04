# Phase 1 findings: the negative set, and why the positives were the real problem

Status as of 2026-08-04. This is the live record of what Phase 1 actually found, as
opposed to what PLAN_v1.md assumed. Read it before touching the training stage.

## 1. The plan assumed order 10^2 independent positives. There is roughly one.

The nine hand-curated positives (IsPETase, ThermoPETase, FAST-PETase,
IsPETase-W159H/S238F, LCC, LCC-ICCG, TfCut2, Cut190, and MHETase separately) collapse
into a **single cluster at both 30% and 50% identity**.

By spec section 8's own rule (splits by sequence cluster, never by sequence) that is one
independent example. No cluster-split evaluation is possible over them, and any
cross-validation reports a number produced entirely by leakage: the engineered variants
are 97.6% to 99.3% identical to their parents, so a variant of the training example sits
in the test fold every time.

This was not a subtle failure. The first ungrouped trivial baseline returned AUC 0.9996.

**Resolution taken:** the ESTHER `Polyesterase-lipase-cutinase` family was harvested from
UniProt in full (83 members, 79 new), taking the positive set to 87 sequences in **11
clusters at 30% identity**, which makes a cluster-split evaluation possible at all.

**The cost, which must be reported in the Methods tab:** those 79 carry ESTHER family
annotation, not measured PET activity. They are stored with `source_ref='ESTHER-family'`
and a note saying so. Spec section 8's last bullet already requires performance to be
reported separately for the quantitatively characterised and annotation-only subsets;
this is now the dominant reason why.

## 2. Risk 1 fired for real, and the diagnosis was not what the spec anticipated

With cluster-grouped cross-validation, the trivial baseline (20 amino-acid fractions plus
length, no fold information, no embedding) scored:

| Negative set | AUC (cluster-grouped) | Null | Verdict |
|---|---|---|---|
| First attempt, families matched on length/identity/taxonomy | 0.954 +/- 0.016 | 0.495 | FAIL |
| Plus secretion matching | 0.845 +/- 0.121 | 0.432 | MARGINAL |

The coefficients diagnosed the first failure precisely: negatives were Leu-rich,
positives Ser/Thr/Gly/Pro-rich. That is a **secreted-versus-cytoplasmic** signature, not
polyester chemistry. Every characterised polyesterase is a secreted enzyme (IsPETase
carries a signal peptide at 1-27, LCC at 1-34, TfCut2 at 1-40), while the allowlisted
negative families (`6_AlphaBeta_hydrolase`, `Epoxide_hydrolase`, `Proline_iminopeptidase`,
`Homoserine_transacetylase` and similar) are largely intracellular.

A model trained against that set could have scored well by learning "is this exported?"
and nothing about PET at all.

**Resolution taken:** negatives are now additionally required to carry a signal peptide
(`ft_signal:*`), a fourth matching axis beyond the three the spec lists. That is what
moved 0.954 to 0.845.

## 3. What is still wrong

0.845 is above the 0.75 pass threshold. Residual separation is still compositional
(P/S/T/G up, K/L down), and the likely remaining confound is **taxonomic GC content**:
the positives are dominated by high-GC Actinobacteria (Thermobifida, Streptomyces,
Amycolatopsis), whose proteins are systematically Ala/Gly/Pro-rich, while the negatives
still include Proteobacteria (Paraburkholderia, Bradyrhizobium).

The spread is also wide (std 0.121, p10 0.71, p90 0.95) because 11 positive clusters is
still very few, so the point estimate should not be over-read in either direction.

### What was then tried, and what it delivered

All three shortlisted options were implemented on 2026-08-04.

| Action | Baseline AUC | Verdict |
|---|---|---|
| Length + identity + genus matching only | 0.954 | FAIL |
| Plus signal-peptide (secretion) matching | 0.845 | MARGINAL |
| Plus phylum matching | **0.842** | MARGINAL |

**Phylum matching contributed essentially nothing (0.845 to 0.842), and that is the
useful result.** The residual separation is therefore *not* driven by GC content, which
was the leading hypothesis. The surviving signal (S/P/T/C up, L/K down) looks like the
intrinsic composition of the polyesterase fold itself: a Ser/Thr-rich surface, cysteines
forming the stabilising disulfides, prolines in the loops.

That reframes the number. At 11 positive clusters with a spread of p10 0.72 to p90 0.96,
the honest reading is not "the negative set is still broken" but "amino-acid composition
carries real, non-trivial information about whether something is a polyesterase". A
composition classifier scoring 0.84 is partly measuring biology, not purely exploiting an
artefact.

**Decisions taken as a result:**

1. The composition baseline is now a **permanently reported metric**, stored in
   `training_runs.composition_baseline_auc` alongside `retrieval_baseline_auc`. It is no
   longer only a pre-training gate. Any claim the learned model makes must clear both
   baselines, which is a stricter and more honest bar than the spec originally set.
2. `training_runs.n_positive_clusters` records independent units rather than the raw
   positive count, so no future run can quietly report "87 positives" when it has 11.
3. Positives now carry an **evidence level** derived from UniProt `protein_existence`:
   23 at `1: Evidence at protein level`, 56 predicted or inferred. Stored as
   `source_ref='ESTHER-family-protein-evidence'` versus `'-predicted'`. Spec section 8's
   requirement to report the characterised and annotation-only subsets separately is now
   mechanically supported rather than aspirational.

Full PAZy curation, with measured rates landing in `activity_measurements`, remains
outstanding and is still the highest-value unblocked work.

### Options, in the order I would try them

1. **Match negatives on phylum, or directly on GC content**, so the composition confound
   is removed rather than reduced. Cheapest, and directly targets the diagnosed cause.
2. **Restrict negatives to secreted actinobacterial esterases only.** Smaller and harder
   set; risks too few negatives to matter.
3. **Accept MARGINAL and proceed**, recording the residual as a known bound on what any
   score can mean. Defensible only if the Methods tab states it plainly and the retrieval
   baseline comparison (spec section 8) becomes the primary claim rather than the AUC.
4. **Invest in real PAZy curation** to get more genuinely characterised, measured
   positives. Highest value, slowest, and the only option that fixes the root cause
   rather than the symptom.

Options 1 and 4 are not alternatives to each other and should both happen.

## 4. Smaller corrections made along the way

- `Q6A0I4` was initially curated as Cut190; UniProt returns **TfCut2** (*Thermobifida
  fusca*). Cut190 is `W0TJ64`, and its strain assignment (AHK190 versus type strain P101,
  both 304 aa) is **still unresolved**.
- MHETase's ESTHER family is `Tannase`, confirming spec section 2 point 5 directly: it
  shares no family with the PETases, so a PETase-seeded profile search cannot reach it.
  It is excluded from this negative-set work and needs its own v2 pipeline.
- A separate `Cutinase` family exists alongside `Polyesterase-lipase-cutinase`. Its 111
  members are stored as near-misses, which is what spec section 5.2 asks for.
- PHB depolymerases and the `Tannase` family are excluded from negatives entirely rather
  than labelled inactive: calling a polyester-active enzyme "no polyester activity" would
  teach the model the opposite of what is wanted.

---

# Phase 2a findings: the triad filter works, and it is profile-bound

## 5. The mechanism validates, because the reference had a known answer

The triad filter aligns candidates to a profile HMM with `hmmalign` and reads the columns
corresponding to IsPETase's own verified S160/D206/H237. Running it on IsPETase itself is
therefore a test with an answer known in advance, and it failed the first time:

```
got      Ser@171=A   Asp@217=P   His@251=G
expected Ser@160=S   Asp@206=D   His@237=H
```

Cause: `hmmalign` writes residues that are insertions relative to the model in **lower
case**, and match states in upper case. Both are real residues and both consume a position
in the sequence's own numbering, but the mapping counted only upper-case characters, so it
walked off by the number of insertions before the target.

**This failure is silent in production.** A broken mapping does not raise: the triad simply
reads as incomplete and the candidate is quietly discarded, so the bug looks exactly like
a strict filter working correctly. Nothing except a reference with a pre-known answer
would have caught it. Now covered by `tests/test_triad.py`.

After the fix:

| Set | Complete triad |
|---|---|
| Positives | 79/87 (91%) |
| Hard negatives | 20/252 (8%) |
| Near misses | **0/111 (0%)** |

Named enzymes all resolve sensibly: IsPETase S160/D206/H237, LCC S165/D210/H242,
TfCut2 S170/D216/H248, Cut190 S176/D222/H254.

## 6. The 0% on near misses is a design problem, not a result

Classic cutinases certainly have catalytic triads. Scoring 0/111 means they are not being
found to lack a triad, they are failing to align to a Polyesterase-lipase-cutinase
profile well enough for the triad columns to map at all.

That matters because the near misses exist precisely to define the decision boundary
(spec section 5.2). A filter that discards them wholesale removes the examples the model
most needs to see, and it does so before scoring, where the discard is invisible.

**The single pooled profile used for this validation is the cause.** `profiles.py` already
implements `build_from_clusters` for one profile per cluster; the validation short-cut to
a single profile over all 87 positives to test the mechanism. Production needs the
per-family design: build a profile per cluster (including a Cutinase profile), align each
candidate to its best-scoring profile, and read the triad from that profile's own
reference rather than forcing everything through IsPETase's numbering.

Until that lands, the reported triad-completeness numbers are conditional on the profile
used, and the 8% figure for hard negatives is a floor rather than an estimate.

## 7. MGnify has very little plastic-associated assembly data

Surveying MGnify for the collection choice returned 60 studies across plastic,
plastisphere, landfill and compost search terms. Filtering to those with **assembly**
analyses (the only kind carrying predicted proteins) leaves **24**, and every one of them
has a single assembly. The largest plastic-relevant studies by sample count are amplicon:
`MGYS00001767` "Plastisphere Targeted Locus (Loci)" has 357 samples and no protein
sequences whatsoever.

The assembly studies that do exist are dominated by compost enrichment cultures rather
than plastisphere or landfill.

This materially changes decision 1 in PLAN_v1.md, which assumed one MGnify study would
supply order 10^6 protein sequences. Options:

1. **Take the compost assemblies anyway.** Compost is a defensible polyesterase habitat
   (LCC itself is leaf-branch compost derived), and the data is immediately available.
   Smallest effort, weakest link to the therapeutic framing.
2. **Go to JGI IMG/M or the marine plastisphere assemblies directly**, outside MGnify.
   More work, better matched to the spec's environment priorities.
3. **Assemble from raw reads.** Correct, and far beyond a v1 compute budget on one Mac.

Recommendation: 1 for a working v1 end to end, with 2 as the first expansion. The recall
stage does not care where the FASTA came from, so switching later costs nothing already
built.

---

# Phase 2b: the first real recall run

Data acquired 2026-08-04: **2,220,462 predicted proteins, 858 MB**, across landfill
(Riverton City dump), marine plastisphere (PRJNA777294) and compost (cattle manure,
ZCTH02). Both shortlisted options were satisfied without JGI IMG/M, because MGnify's TPA
assemblies carry the marine plastisphere data and IMG/M would have needed manual
credentials.

## 8. The funnel

| Stage | Surviving | Note |
|---|---|---|
| Scanned | 2,220,462 | |
| MMseqs2 prefilter, E <= 1e-5 | 1,291 | 0.06% |
| Matched a profile (hmmscan) | 638 | |
| Complete catalytic triad | 114 | 18% of profile-matched |
| Unique candidates written | **110** | content-addressed, so 4 were the same protein recovered from two assemblies |

**0.005% of input is retained.** That is aggressive, and it is worth being explicit that
this is a *choice* (E <= 1e-5 against only 87 positives, then a strict triad requirement)
rather than a property of the data. It lands almost exactly where PLAN_v1.md wanted for
v1: 110 candidates, of which the top 50 get structures.

By environment: compost 58, marine plastisphere 40, landfill 12.

Runtime was 1,902 s (32 minutes) for the whole 2.2M on the M1 Max, against the plan's
estimate of 2 to 4 hours for MMseqs2 alone. The estimate was pessimistic.

## 9. The top hits are rediscoveries, and that is the point

The highest-scoring candidates sit at 96 to 99% identity to already-characterised enzymes
(the best is 96% identical to TfCut2). Those are not discoveries, they are the pipeline
proving it can find what is already known.

This is spec section 2's thesis showing up in the data on the first run: **a homology
search ranks the well-known enzymes to the top.** Mean identity to the nearest
characterised enzyme across all 110 candidates is 0.398, so the genuinely interesting
candidates are the low-identity ones that E-value rank pushes *down*. Re-ranking those is
what the learned model is for, and the retrieval numbers now stored on every candidate
(`recall_evalue`, `recall_bitscore`, `recall_profile_identity`) are the baseline it has
to beat.

## 10. A metadata bug worth recording

The first run inferred `source_environment` from filename prefixes and left 36 of 110
candidates labelled `unknown`. Nothing failed: the rows were written, the counts looked
fine, and the error would only have surfaced as an empty environment filter in the
Catalogue and a mis-coloured Home hero plot.

Environment is a ranking-relevant axis in this project, not decoration, so it is now
resolved from each study's MGnify biome via `mgnify.STUDY_BIOME` rather than guessed from
a filename.

---

# Phase 1.1 revisited: activity curation without typing rates out of PDFs

PAZy has no API, so the obvious route to measured activity is transcribing numbers from
papers. That is exactly where fabrication risk lives, and a wrong rate is undetectable
downstream. UniProt carries the same information in machine-readable, citable form, so
every row below is sourced and checkable.

## 11. EC 3.1.1.101 is the label that was missing

`EC 3.1.1.101` is **poly(ethylene terephthalate) hydrolase**: a curator's assignment of
PET-hydrolysing function, not a family guess. UniProt holds 459 entries carrying it
(EC 3.1.1.102 for MHETase holds 62).

**But the tier matters, and nearly got misreported.** Only the 10 Swiss-Prot reviewed
entries carry `ECO:0000269` (experimental evidence from a publication). The other 449 are
TrEMBL entries whose EC comes from `ECO:0000256`, automatic annotation by similarity.
Those were initially labelled `EC-3.1.1.101-unreviewed`, which reads as a mere curation
backlog rather than the substantive difference it is. They are now `EC-auto-annotated`,
with a note saying plainly that the EC number is a prediction of PET activity and not a
measurement.

Positive set after curation: **529, of which 16 are experimentally evidenced.**

| Tier | n | What it means |
|---|---|---|
| `EC-auto-annotated` | 449 | EC 3.1.1.101 by similarity (ECO:0000256). A prediction |
| `ESTHER-family-predicted` | 50 | Family membership only |
| `ESTHER-family-protein-evidence` | 14 | Family, protein observed |
| `EC-experimental` | 10 | EC 3.1.1.101 with ECO:0000269 and PubMed citations |
| Curated wild types and variants | 6 | Hand-curated, sequence-verified |

## 12. 47 measurements, every one carrying its PubMed IDs

| Parameter | Rows | With a parsed value |
|---|---|---|
| Km | 21 | 21 |
| Catalytic activity (PET) | 10 | n/a, qualitative |
| Temperature optimum | 8 | 8 |
| pH optimum | 8 | 8 |

Optima are stored with their prose intact in `raw_text`, because the number alone loses
the substrate and conditions it depends on ("Optimum pH is 8.5 with pNP-butyrate as
substrate"). `comparable_group_id` keys on parameter plus substrate, so a Km on
pNP-butanoate can never be pooled with one on PET film. UniProt often embeds the
conditions in the substrate string, which fragments the groups further, and that is
correct: `km:pnp-butanoate (at 50 degrees celsius and ph 8)` genuinely is not comparable
with `km:pnp-butanoate (at 25 degrees celsius and ph 7)`.

## 13. The therapeutic gap, now measured rather than asserted

Spec section 1.1 tabulates the industrial-versus-therapeutic mismatch as a premise. The
extracted optima turn it into data:

| Enzyme | Topt | pH opt |
|---|---|---|
| Q47RJ6, Q47RJ7 (*T. fusca*) | 60 °C | 8.0 |
| G8GER6 (TfCut1), TfCut2 | 55 °C | 8.0 |
| LCC | 50 °C | 8.5 |
| F7IX06, D4Q9N1 (*T. alba*) | 50 °C | 6.0 |
| **IsPETase** | **40 °C** | 9.0 |

Every characterised PET hydrolase with a measured optimum sits at 40 to 60 °C and mostly
alkaline. IsPETase at 40 °C is the closest thing to a physiological enzyme in the entire
characterised set, and it is still 3 °C above body temperature at a pH optimum of 9.

This is also real training data for the v2 Topt head, which previously had none.

## 14. The gate moved

Re-running the risk-1 trivial baseline with the curated positive set:

| Positive set | Clusters | Composition baseline AUC | Verdict |
|---|---|---|---|
| Curated only | 1 | 0.9996 (leakage) | invalid |
| Plus ESTHER family | 11 | 0.842 | MARGINAL |
| **Plus EC 3.1.1.101** | **29** | **0.778** | MARGINAL |

0.778 against a 0.75 pass threshold, with p10 at 0.663. Still not a pass, but the trend
confirms the Phase 1 diagnosis: a large part of the apparent compositional shortcut was a
small-sample, low-diversity artefact, and it shrinks as real diversity is added. Twenty-nine
independent clusters is also the first point at which the evaluation protocol in spec
section 8 has enough units to mean much.

**Caveat to carry forward:** the EC-annotated set widens the length range to 63 to 835 aa,
which includes fragments and probable multi-domain proteins. That should be filtered before
training, and it is why the matched negative set shrank to 131.

---

# Phase 4: filtering and embedding

## 15. Recall re-run against the curated positives

The 110 candidates from the first recall run were destroyed by a database rebuild during
the curation work (my error: the FASTA files live outside the repo so only compute was
lost). Re-running turned out to be the better path, because the profile library could then
be built from the curated 529 positives across 29 clusters rather than the original 87
across 11.

| | First run | Re-run |
|---|---|---|
| Library built from | 87 positives, 11 clusters | 529 positives, 29 clusters |
| Candidates | 110 | **128** |
| Runtime | 1,902 s | 1,424 s |
| Environments | 36 labelled `unknown` | all resolved |

By environment: compost 69, marine plastisphere 44, landfill 15. The wastewater assembly
(ERZ795023, 26,631 proteins) yielded nothing, which is a reasonable null: it is the only
non-plastic-associated, non-compost source in the set.

## 16. Filtering: marked, not deleted

| Reason | n |
|---|---|
| Included | 720 |
| Length outside 200 to 450 aa | 57 |
| UniProt `Fragment` flag | 8 |
| No sequence (mutation set unconfirmed) | 5 |

Fragments come from UniProt's own flag rather than a length cutoff, because a short
sequence is not necessarily a fragment and a fragment is not necessarily short. The length
window is derived from the experimentally evidenced positives (262 to 319 aa excluding
MHETase), not chosen for roundness.

Exclusions are recorded on the row. The catalogue is a deliverable in its own right, so a
fragment still belongs in it: what it must not do is silently become a training example.

## 17. Embedding, and one number worth carrying into training

848 vectors at 480 dimensions: 720 training sequences (30 s) plus 128 candidates (9 s),
nothing truncated. The plan budgeted 25 to 40 minutes for ~10^4 sequences; at ~20 seq/s on
CPU that estimate holds.

Pooling excludes padding as well as CLS and EOS. Pooling over padding would make a vector
depend on which other sequences shared its batch, and that class of bug is dangerous
precisely because the resulting numbers still look reasonable.

Sanity checks pass, and one of them is a warning about what comes next:

- FAST-PETase to IsPETase cosine **1.000** (they differ by 5 residues in 290)
- LCC to LCC-ICCG **0.999**
- All characterised PET enzymes sit at **0.96 or above** to each other
- Candidate to nearest known enzyme: max 1.000, **median 0.931**, min 0.771

Everything of interest lives in a very tight region of the embedding space. The head is
therefore discriminating small differences inside a dense cluster, not separating
well-spaced groups, which is a further reason to trust calibration curves over headline
AUC (spec section 8).

---

# Phase 5: the labelling experiment, and why the head cannot yet be trusted

Three labelling schemes, cluster-grouped at 30% identity, on frozen ESM-2 embeddings.

| Scheme | Positives | Unlabelled | Negatives | Clusters with positives | Folds | AUC |
|---|---|---|---|---|---|---|
| naive (annotation as truth) | 500 | 0 | 220 | 25 | 4 | **1.000 +/- 0.000** |
| evidence only (ECO:0000269) | 13 | 0 | 220 | **1** | **0** | not evaluable |
| PU (Elkan-Noto) | 13 | 487 | 220 | **1** | **0** | not evaluable |

## 18. AUC 1.000 is the failure, not the success

The naive head separates the 500 annotated positives from the 220 hard negatives
perfectly, across every cluster-grouped fold, and beats the composition baseline decisively
(1.000 against 0.778).

That is not evidence it learned PET chemistry. **449 of those 500 labels were themselves
assigned by sequence similarity** (`ECO:0000256`), so a sequence-embedding model
reproducing them is close to tautological: the head recovered the annotation rule, which
was a similarity rule, using a similarity-based representation. The task as labelled is
circular, and a perfect score is what circularity looks like.

This is exactly the failure spec section 5.3 warns about from the other direction, and it
would have been invisible without deliberately splitting the labels by evidence.

## 19. The evidenced positives are one cluster, again

Only 13 of the filtered training positives carry experimental evidence, and **all 13 fall
into a single cluster at 30% identity**. StratifiedGroupKFold therefore cannot construct a
fold with evidenced positives in both train and test, so both the evidence-only and PU
schemes produce zero valid folds.

This is the Phase 1 finding recurring at a different scale. It is not an artefact of the
filtering: the experimentally characterised PET hydrolase field genuinely is one tight
sequence neighbourhood (IsPETase, LCC, TfCut1/2, Cut190, *T. alba* and *T. cellulosilytica*
cutinases all sit together).

**Consequence: a head trained only on measured PET activity cannot currently be evaluated
by the protocol spec section 8 requires.** Not "scores badly": cannot be scored at all.

## 20. The two schemes disagree exactly where it matters

Despite the evaluation problem, both heads can still be fitted and used to rank the 128
candidates:

| | naive | evidence only |
|---|---|---|
| Mean score | 0.564 | 0.234 |
| Median score | 0.519 | 0.094 |
| Candidates scored above 0.5 | 71 / 128 | 26 / 128 |
| Spearman rho between them | 0.882 | |
| **Top-10 overlap** | **3 / 10** | |

Broad agreement (rho 0.88) with sharp disagreement at the top, which is precisely the part
that matters: the structure budget takes the top 50, and the two schemes would send
substantially different candidates for folding.

Both also correlate with the retrieval baseline they are supposed to improve on
(Spearman rho 0.569 naive, 0.622 evidence-only against -log10 E-value), so a large part of
what either head is doing is reproducing the ranking retrieval already gave us.

## 21. What this settles

The experiment was run to decide whether extended curation or structural features should
come next. It answered clearly: **labels are the bottleneck, and the specific thing needed
is evidenced positives from sequence clusters other than the one the field already knows.**

That is a hard requirement, not a matter of effort. Options, honestly ranked:

1. **Extend the evidenced set outward.** PAZy and the wider plastic-degrading literature
   cover polyesterases beyond the PET-cutinase cluster (PHB depolymerases, polyurethanases,
   polyamidases, MHETase). Broadening the target from "PET hydrolase" to "polyesterase"
   would buy cluster diversity at the cost of specificity, which is a real trade and should
   be a deliberate decision rather than a drift.
2. **Report the ranking honestly as retrieval plus a family prior**, and drop the claim of
   learned PET discrimination for v1. Spec section 2 already anticipates this: "if the
   model cannot beat E-value rank on held-out characterised enzymes, that is a cheap and
   publishable negative result." That is where the evidence currently points.
3. **Use the structure stage as the discriminator instead.** Cleft width and the aromatic
   clamp are the features the spec expects to separate real polyesterases from
   soluble-ester esterases, and they are measured from geometry rather than inherited from
   an annotation. This route does not depend on more labels.

Option 3 is the only one that adds genuinely independent signal, and it is unblocked: the
128 candidates can be folded and measured without any further curation.

---

# Phase 6: geometry, the first signal that does not inherit an annotation

## 22. Risk 7 cleared: the geometry recovers published triads

Validated against experimental structures before touching any prediction. The triad is
found by **geometry**, not by sequence position: the code searches for a spatially
connected Ser/His/Asp, so three residues present in sequence but not in contact are
correctly not a triad.

| PDB | Enzyme | Triad found | Expected | Ser-His | His-Asp | Connected |
|---|---|---|---|---|---|---|
| 6EQE | IsPETase | S160/D206/H237 | S160/D206/H237 | 2.94 Å | 3.07 Å | yes |
| 4EB0 | LCC | S165/D210/H242 | S165/D210/H242 | 2.45 Å | 2.75 Å | yes |
| 4CG1 | TfCut2 | S130/D176/H208 | S170/D216/H248 | 2.77 Å | 3.03 Å | yes |
| 1CEX | *F. solani* cutinase | S120/D175/H188 | not annotated | 2.71 Å | 3.19 Å | yes |

The 4CG1 row looks like a failure and is not: the offset is **exactly 40 on all three
residues**, which is precisely TfCut2's signal peptide length (1 to 40). The PDB entry
numbers the mature protein while UniProt numbers the precursor. A constant offset across
all three is the signature of a numbering convention rather than a detection error, and
`tests/test_geometry.py` pins it so nobody "fixes" it later.

## 23. Cleft width separates polyesterases from cutinases, as the spec predicted

| PDB | Enzyme | Cleft width | Depth | Aromatic clamp |
|---|---|---|---|---|
| 6EQE | IsPETase | **20.90 Å** | 2.56 | PHE201, TRP159, **TRP185**, TYR87 |
| 4EB0 | LCC | 17.82 Å | 2.20 | PHE125, PHE243, TRP190, TYR127, TYR95 |
| 4CG1 | TfCut2 | 12.77 Å | 2.29 | PHE209, TRP155, TYR60 |
| 1CEX | *F. solani* cutinase | **11.17 Å** | 1.29 | PHE147, TYR119, TYR149, TYR191 |

Spec section 6: "Cleft width and the aromatic clamp are the features most likely to
separate real polyesterases from soluble-ester-only esterases, so surface them
prominently." That prediction holds on the first test. The range is nearly two-fold and it
orders correctly: IsPETase widest, the classic fungal cutinase narrowest.

The clamp IsPETase returns includes **TRP185**, the mobile tryptophan whose flexibility
spec section 2 point 1 names as the distinguishing feature against cutinases with a fixed
equivalent, and TYR87. Neither was supplied to the code: both fall out of a distance
search around the nucleophile.

## 24. Why this matters more than its place in the plan suggests

Every sequence-derived label in this project traces back to somebody else's annotation,
and Phase 5 showed what that costs: a head trained on those labels scored AUC 1.000 by
reproducing a similarity rule with a similarity representation.

Cleft width is different in kind. It is measured off coordinates, so it is the first
feature here that is independent of the annotation, and on the four structures where the
answer is known it already orders them by the property of interest. That is the route out
of the circularity, and it needs no further curation.

---

# Phase 6 continued: ESMFold works, and section 23 above was partly wrong

## 25. CORRECTION to section 23

Section 23 reported cleft widths of IsPETase 20.90, LCC 17.82, TfCut2 12.77 and
*F. solani* cutinase 11.17 Å, and called it "nearly a two-fold range, ordered correctly".
**That ordering was an artefact of a brittle metric, and the claim was overstated.**

The measure took a MAX over pairwise aromatic separations within a hard 12 Å radius of the
nucleophile. Folding IsPETase with ESMFold and measuring its own prediction exposed it:

| | Cleft width |
|---|---|
| 6EQE crystal | 20.90 Å |
| ESMFold prediction of the same protein | **11.09 Å** |

The fold was fine. PHE201 sits 11.7 Å from the nucleophile in the crystal and 12.1 Å in
the prediction, so a **0.4 Å coordinate shift moved one residue across the cutoff**, and
because the width was a max over pairs, that single membership change halved the answer.
It would have placed the reference PETase in cutinase territory. Both numbers looked
perfectly plausible, which is what makes this class of bug dangerous.

Fixed by widening the radius to 16 Å and taking the **90th percentile** of pairwise
separations rather than the max, so no single residue can decide the value.
Crystal-to-prediction disagreement fell from 9.81 Å to 1.81 Å.

**What survives, stated honestly:**

| Enzyme | Cleft width (p90) |
|---|---|
| LCC | 24.22 Å |
| TfCut2 | 22.43 Å |
| IsPETase | 21.27 Å |
| *F. solani* cutinase | **14.93 Å** |

The classic fungal cutinase is clearly narrower than the Polyesterase-lipase-cutinase
family members, by around 7 Å. But those members **do not order among themselves by PET
activity**: IsPETase, the best ambient-temperature PET degrader of the three, has the
narrowest cleft of them.

So the demonstrated result is a **family-level separation**, not a PET-activity one.
Separating a fungal cutinase from the PLC family is a much easier task than separating
PET-active from PET-inactive within the PLC family, and only the former is shown. That
weakens, but does not eliminate, the case for geometry as the route out of the Phase 5
circularity: it is still an annotation-independent measurement, it just has not yet been
shown to track the property of interest.

## 26. ESMFold: risk 4 cleared, with a caveat about size

| | |
|---|---|
| Model | `facebook/esmfold_v1`, a single **8.44 GB** weights file (it bundles ESM-2 **3B**) |
| Download | ~50 min on this connection |
| Load time | 1,849 s (31 min), 11.3 GB resident |
| Fold time | **109 s for 290 residues** (0.38 s/residue) |
| Extrapolated for 50 candidates | **~94 minutes** |
| Mean pLDDT | 0.90 on a 0 to 1 scale (equivalent to 90) |

Folding is fast enough; the cost is loading. No fallback to Boltz-2 is needed.

**The prediction reproduces the crystal active site exactly**: S160/D206/H237, connected,
Ser-His 2.98 Å and His-Asp 2.78 Å against the crystal's 2.94 and 3.07 Å. The aromatic
clamp comes back with TRP185 and TYR87.

Note for anyone reading the code: ESMFold's pLDDT is on a **0 to 1** scale here, not
0 to 100. The smoke test first printed "mean pLDDT 1.0", which looked like a broken tensor
read and was actually a rounded CA-atom mean on the unfamiliar scale.
