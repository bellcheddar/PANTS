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
