# PANTS v1 Implementation Plan

**Planning pass output for `pants_spec_v1.md` section 13.**
**Date: 2026-08-04. Target deployment: `pants.mdeller.com` (port 8005).**

---

## 0. What the codebase actually says (corrections to the spec)

Three places where the spec's assumptions do not match evidence gathered from AlphaFraud, BoltzMaker, FlexAppeal, chem_sage-web, chatPDB-web and TopPDBLX.

**Section 3.1's "under about 2 GB" reads as if PANTS is the only tenant.** It is not. The droplet already runs five gunicorn services: AlphaFraud (port 8000, 3 sync workers), chem_sage-web (8001, 1 gevent worker, thin proxy to a remote HF Space), chatPDB-web (8002, same pattern), BoltzMaker-web (8003, 3 sync workers, thin flask-only venv), FlexAppeal (8004, 2 sync workers, `MemoryMax=1200M` cgroup cap because MDTraj can spike). FlexAppeal's own `PROJECT_PLAN.md` states it plainly: "the droplet has 3.8 GB of RAM shared with four other applications." PANTS is app number six, on port 8005. The real question is not "does 2 GB fit in 3.8 GB", it is "does PANTS fit in whatever is left after five other apps and the OS", which is a materially smaller number (worked in section 5).

**Section 7's Home hero plot (activity probability versus predicted Topt) assumes a trained Topt head.** Section 13's own v1 scope line promises only "one calibrated activity head" (PET activity). There is no Topt head in v1. Shipping a fake or placeholder Topt number to make the spec'd chart work would be dishonest by the app's own Methods-tab standard. Resolution in section 7.

**Local disk is tighter than "an M1 Max" suggests.** `df -h` on the Data volume shows 3.6 TiB total, **62 GiB free** (verified 2026-08-04). That is the real local budget. The sibling project `TopPDBLX/PHASE0_PLAN.md` already hit this and routes around it: `data/raw` and `data/interim` are symlinks to `~/TopPDBLXData/`, outside iCloud sync. PANTS must do the same, and decision 1 (metagenome collection) has to be sized against 62 GB, not a cluster-scale assumption.

**Verified on this machine, 2026-08-04:** mmseqs2 present at `/opt/homebrew/bin/mmseqs`. HMMER **not yet installed** (`hmmsearch`/`hmmbuild` absent); the arm64 bottle exists, so Phase 0.3's `brew install hmmer` is expected clean but is a real step, not a no-op. System Python is 3.14.3.

---

## 1. Task breakdown, dependency order, critical path

Phased, each task sized for a solo developer using Claude Code, working days unless marked (h) for hours. **Critical path in bold.**

### Phase 0: Scaffold (blocks everything, ~1.5 days)

- **0.1 Repo layout** (2h): create `pipeline/`, `app/`, `models/`, `data/`, `eval/`, `deploy/`, `tests/` directly under `PANTS/` (repo root = namespace, exactly as AlphaFraud's repo root holds `AlphaFraud.py` and `alphafraud/` side by side, and exactly as spec section 11 shows `pants/` as the literal root, not a nested package). Add `PANTS.py` CLI entry mirroring `AlphaFraud.py`'s `init`/`run`/`serve` command shape.
- **0.2 Two venvs, not one** (2h, ports FlexAppeal/BoltzMaker-web's proven split): a **web venv** (flask, gunicorn, jinja2, python-dotenv, gemmi: no torch, ever, in v1) and a **pipeline venv** (torch, transformers, scikit-learn, biotite, gemmi, tmtools, pandas). Smoke-test `pip install torch` in a scratch Python 3.14 venv first (30 min), since 3.14.3 is this house's proven interpreter for AlphaFraud/BoltzMaker-web but chem_sage deliberately pinned 3.11+ for its own torch stack. That discrepancy is unexplained and needs a direct check, not an assumption. Fallback: pipeline venv on 3.11/3.12 if 3.14 wheels are not there yet.
- **0.3 Install HMMER, confirm MMseqs2** (15 min): `brew install hmmer` (arm64 bottle, 3.4). mmseqs2 already present (18-8cc5c).
- **0.4 Data dir outside iCloud** (30 min): `mkdir ~/PANTSData`, symlink `data/raw` and `data/interim` into it, exactly as TopPDBLX does for the same reason.
- **0.5 `db.py` and `manifest.py`** (0.5 day): port AlphaFraud's `connect()`/WAL/`_retry_write` pattern and TopPDBLX's manifest convention (every stage writes `manifests/<stage>_<timestamp>.json`: input SHA256s, output SHA256s, tool versions, git commit, schema version, wall time). Full schema in section 3.

### Phase 1: Data acquisition (depends on 0)

- **1.1 PAZy positive curation** (1 day, timeboxed): the named enzymes in spec section 1 plus whatever PAZy's own list yields. PAZy has no documented API: this is manual/LLM-assisted curation, not a scrape.
- **1.2 ESTHER hard negative harvest** (1.5 to 2 days): family list, then sequence pull via UniProt family queries, reusing AlphaFraud's `http.py` (requests + tenacity, already proven in production).
- **1.3 Characterised reference set** (parallel with 1.2, 1 day): FAST-PETase, IsPETase, DuraPETase, HotPETase, LCC-ICCG, Z1-PETase, TurboPETase, Cut190\*\*SS, HGMP01, MHETase: sequences, PDB IDs, activity notes.
- **1.4 Metagenome collection selection and download** (0.5 to 1 day, size-checked first, see risk 5): **one** MGnify plastisphere/landfill study, not breadth-first (arithmetic in section 4).
- **1.5 PDB reference structures** (2h): IsPETase 6EQE and relatives, via AlphaFraud's `pdb.py` (RCSB GraphQL, proven).

### Phase 2: Recall (depends on 1.1 to 1.4) — **critical path**

- **2.1 Build profile HMMs** (0.5 day): `hmmbuild` from PAZy positives plus ESTHER families.
- **2.2 MMseqs2 prefilter** (0.5 day dev plus compute, section 4) over the metagenome protein set.
- **2.3 HMMER sensitive pass** (0.5 day dev plus compute) on MMseqs2 survivors.
- **2.4 Triad/oxyanion-hole completeness filter** (0.5 day): write survivors into `candidates`, report discard counts (spec section 5.1 requires this).

### Phase 3: Negatives and PU setup (overlaps Phase 2, depends on 1.2)

- **3.1 Length/identity/taxonomy matching** (0.5 to 1 day).
- **3.2 Near-miss curation** (0.5 day, manual: cutinases active on soluble esters, not PET).
- **3.3 PU class prior decision** (0.5 day: v1-blocking, see section 6).

### Phase 4: Embedding (depends on 2, 3) — **critical path**

- **4.1 ESM-2 t12-35M batch embed** (0.5 day dev plus ~0.5 to 1.5h compute, section 4). CPU only, not MPS (justification in section 2).

### Phase 5: Model training (depends on 4) — **critical path**

- **5.1 PET activity head** (1 to 1.5 days): logistic regression on frozen embeddings, Elkan-Noto PU correction, `CalibratedClassifierCV` (Platt/isotonic).
- **5.2 Evaluation** (1 day): MMseqs2 cluster splits at 30% and 50%, leave-one-family-out, retrieval baseline versus HMMER E-value rank, reliability diagram, prospective holdout by date.
- **5.3 Write scores** (0.5 day) into `scores`.

### Phase 6: Structure layer (depends on 5 for top-50 ranking) — **critical path**

- **6.1 Select top 50** by activity probability (1h).
- **6.2 ESMFold batch prediction** (dev 0.5 day plus compute, section 4). Smoke-test ONE structure first (risk 4).
- **6.3 Superpose onto 6EQE** (0.5 day): reuse `tmtools` exactly as AlphaFraud's `compare.py` does. **Note the thread-safety gotcha**: TM-align's C bindings use static buffers and are not thread-safe; AlphaFraud guards this with `_TMALIGN_LOCK`. Reuse the lock, or parallelise by process not thread.
- **6.4 Active-site geometry extraction** (1 to 1.5 days, genuinely new code, nothing to port): triad distances and angles, oxyanion hole, cleft width and depth, aromatic clamp, solvent accessibility. Validate against known PDB structures first (risk 7).
- **6.5 mmCIF assets plus `structures`/`geometry` tables** (0.5 day).

### Phase 7: Web app (scaffolding can start as soon as Phase 0 lands, in parallel with Phases 1 to 6, using fixture data)

- **7.1 App factory and brand theme** (0.5 day): port `alphafraud/static/brand.css`, house header.
- **7.2 Home tab** (0.5 day, see section 7 for hero-plot substitute).
- **7.3 Catalogue tab** (1 day): Tabulator table, filters, composite slider (single-axis in v1, scaffolding for v2).
- **7.4 Candidate tab** (1.5 days): Mol\* viewer over static mmCIF, per-residue pLDDT track, score panel, sequence with triad annotated, nearest characterised relatives.
- **7.5 Methods tab** (1 day): mostly templating over `manifests`/`training_runs`, limitations section lifted verbatim from spec section 9.
- **7.6 Mol\* wrapper** (0.5 day, shared component reused across tabs).

### Phase 8: Deploy (depends on 7) — **critical path tail**

- **8.1 Mirror AlphaFraud/FlexAppeal `deploy/`** (0.5 to 1 day): `gunicorn.conf.py` (port 8005), `nginx-pants.conf`, `pants-web.service`, `pants-run.service` and `.timer`, `provision.sh`, `deploy.sh`.
- **8.2 Provision, DNS, certbot** (0.5 day, some idle wait).
- **8.3 `mdeller-landing/apps.json` entry plus `./deploy.sh`** (2h: one JSON block, confirmed against the live file).
- **8.4 `free -m` on the real droplet before going further** (15 min, see risk 6: this replaces every estimate in section 5 with ground truth).
- **8.5 End-to-end smoke test** (0.5 day).

**Critical path total: roughly 3 to 4 working weeks (15 to 20 days)**, dominated by Phase 1's manual curation (biology judgement, not code) and Phase 6.4's novel geometry code. Overlapping Phase 7's UI scaffolding with Phases 1 to 6 is the main schedule lever available.

---

## 2. Concrete Python packages per stage

| Stage | Packages | Apple Silicon note |
|---|---|---|
| Recall | `hmmer` (brew, CLI subprocess), `mmseqs2` (brew, CLI subprocess), `biopython` | mmseqs2 **verified installed** (18-8cc5c); hmmer arm64 bottle 3.4 available, install in 0.3. Shell out via subprocess, matching BoltzMaker's own tool-invocation house style, rather than adding `pyhmmer` bindings. |
| Negatives | `requests`, `tenacity`, `pandas`, `mmseqs2` (identity matching) | Reuse AlphaFraud's `http.py` verbatim (proven in production since 2026-07-14). |
| Embed | `torch>=2.x` (CPU only), `transformers>=4.46` (same pin as chem_sage) | **CPU, not MPS**, deliberately: ESM-2 t12-35M is 35M params, cheap enough that MPS's dtype quirks and documented silent-CPU-fallback ops (BoltzMaker hit this with `torch.linalg.svd`) are not worth the complexity. Save MPS effort for the structure stage where it matters. |
| Train | `scikit-learn>=1.4` (`CalibratedClassifierCV`), `numpy`, `scipy` | No PU-learning package needed: Elkan-Noto correction is ~50 lines of numpy. This matches an existing house convention (AlphaFraud's `db.py` hand-rolls a 2-parameter Newton-Raphson logistic regression in `cw_rate_trend()` rather than pulling in `statsmodels`). No W&B: it is not used anywhere in this stack outside MLX fine-tuning. Log runs to the `training_runs` table instead, mirroring AlphaFraud's `runs` table. |
| Structure | `transformers` (`facebook/esmfold_v1`) primary; BoltzMaker subprocess fallback; `gemmi` (mmCIF I/O, FlexAppeal's choice); `biotite` (alignment/torsion, AlphaFraud's choice); `tmtools` (TM-align, reuse directly) | **Needs real verification, not assumption.** HF's `EsmForProteinFolding` reimplements the layers the original ESMFold repo needed the CUDA-only `openfold` package for, so it should run on CPU/MPS without a compiled-kernel build, but smoke-test on day 1 of Phase 6 (risk 4). BoltzMaker's README documents real MPS pain at scale (hard crash above ~1250 residues from an unchunked matmul, since patched; `torch.linalg.svd` has no MPS kernel and silently falls back to CPU). PETase-like targets are ~250 to 300 residues, monomeric, well clear of that ceiling either way. |
| Harmonise | `pandas` | Manual/LLM-assisted extraction for v1's modest scope (10 to 30 rows); no scraping package needed. |
| App (web venv) | `flask>=3.0`, `gunicorn>=21.2`, `jinja2`, `python-dotenv` | **Identical pins to AlphaFraud's `requirements.txt`, proven.** No torch or transformers in this venv at all: that is the whole point of the two-venv pattern BoltzMaker-web established, so the always-on gunicorn process stays light. |

---

## 3. SQLite schema

Mirrors AlphaFraud's `db.py` conventions: WAL mode set once in `init_schema()`, `connect()`/`_retry_write` wrapper, large JSON blobs kept in sidecar tables so hot list and sort queries stay fast. AlphaFraud's own experience is worth repeating here: an uncovered index on a 28-column table cost 62 ms per page load until a covering index fixed it. Replicate that discipline.

```sql
-- ============================================================
-- Provenance / manifest / versioning (spec section 12 requirement)
-- ============================================================
CREATE TABLE IF NOT EXISTS runs (
    id           INTEGER PRIMARY KEY AUTOINCREMENT,
    stage        TEXT NOT NULL,          -- recall | negatives | embed | train | structure | harmonise
    label        TEXT NOT NULL,          -- run date or dataset version tag
    started_at   TEXT,
    finished_at  TEXT,
    n_input      INTEGER DEFAULT 0,
    n_output     INTEGER DEFAULT 0,
    n_discarded  INTEGER DEFAULT 0,
    status       TEXT DEFAULT 'running', -- running | done | error
    params_json  TEXT
);

CREATE TABLE IF NOT EXISTS manifests (
    id                  INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id              INTEGER,
    stage               TEXT NOT NULL,
    input_hashes_json   TEXT,            -- {path: sha256}
    output_hashes_json  TEXT,
    tool_versions_json  TEXT,            -- {hmmer: '3.4', mmseqs2: '18-8cc5c', transformers: '4.46', ...}
    model_version       TEXT,
    git_commit          TEXT,
    wall_time_s         REAL,
    written_at          TEXT,
    FOREIGN KEY(run_id) REFERENCES runs(id)
);

CREATE TABLE IF NOT EXISTS data_sources (
    name          TEXT PRIMARY KEY,      -- PAZy | ESTHER | PDB | MGnify:<study_id> | Meltome | FireProtDB
    version       TEXT,
    retrieved_at  TEXT,
    n_records     INTEGER,
    license       TEXT,                  -- populated from v1 so decision 7 never needs re-derivation
    source_url    TEXT
);

-- ============================================================
-- Candidates (recall output, catalogue backbone)
-- ============================================================
CREATE TABLE IF NOT EXISTS candidates (
    candidate_id                  TEXT PRIMARY KEY,  -- 'PANTS-' + first 12 hex of SHA256(sequence): stable across reruns
    source_environment            TEXT,              -- plastisphere | landfill | marine | wastewater | compost
    assembly_id                   TEXT,
    contig_id                     TEXT,
    sample_accession              TEXT,
    taxon_name                    TEXT,
    ncbi_taxid                    INTEGER,
    taxonomy_lineage              TEXT,
    sequence                      TEXT NOT NULL,
    seq_length                    INTEGER,
    has_complete_triad            INTEGER,           -- 0/1, filter gate per spec section 5.1
    triad_positions_json          TEXT,
    oxyanion_hole_positions_json  TEXT,
    recall_method                 TEXT,              -- mmseqs2 | hmmer
    recall_evalue                 REAL,
    recall_bitscore               REAL,
    recall_profile_identity       REAL,
    nearest_characterised_id      TEXT,
    nearest_characterised_identity REAL,
    discovered_run_id             INTEGER,
    first_seen_at                 TEXT,
    FOREIGN KEY(discovered_run_id) REFERENCES runs(id),
    FOREIGN KEY(nearest_characterised_id) REFERENCES characterised_enzymes(enzyme_id)
);

CREATE INDEX IF NOT EXISTS idx_candidates_env   ON candidates(source_environment);
CREATE INDEX IF NOT EXISTS idx_candidates_triad ON candidates(has_complete_triad);
CREATE INDEX IF NOT EXISTS idx_candidates_taxon ON candidates(taxon_name);

-- ============================================================
-- Scores (model output; sidecar of candidates)
-- ============================================================
CREATE TABLE IF NOT EXISTS scores (
    candidate_id             TEXT PRIMARY KEY,
    pet_activity_prob        REAL,
    pet_activity_calibrated  INTEGER,   -- 0/1: Platt/isotonic applied
    mhet_activity_prob       REAL,      -- NULL in v1 (v2 head)
    topt_pred_c              REAL,      -- NULL in v1 (v2 head)
    topt_lower               REAL,
    topt_upper               REAL,
    expressibility_prob      REAL,      -- NULL in v1
    aggregation_score        REAL,      -- NULL in v1
    therapeutic_composite    REAL,      -- v1: equals pet_activity_prob (single-axis)
    model_version            TEXT,
    scored_at                TEXT,
    FOREIGN KEY(candidate_id) REFERENCES candidates(candidate_id)
);

CREATE INDEX IF NOT EXISTS idx_scores_composite ON scores(therapeutic_composite);
CREATE INDEX IF NOT EXISTS idx_scores_activity  ON scores(pet_activity_prob);

-- ============================================================
-- Structures and active-site geometry (spec section 6)
-- ============================================================
CREATE TABLE IF NOT EXISTS structures (
    candidate_id            TEXT PRIMARY KEY,
    structure_method        TEXT,       -- esmfold | boltz2
    mmcif_path              TEXT,
    plddt_mean              REAL,
    plddt_per_residue_json  TEXT,
    tm_score_to_ispetase    REAL,
    rmsd_ca_to_ispetase_A   REAL,
    superposition_reference TEXT DEFAULT '6EQE',
    predicted_at            TEXT,
    model_version           TEXT,
    FOREIGN KEY(candidate_id) REFERENCES candidates(candidate_id)
);

CREATE TABLE IF NOT EXISTS geometry (
    candidate_id                  TEXT PRIMARY KEY,
    triad_ser_resnum              INTEGER,
    triad_his_resnum              INTEGER,
    triad_asp_resnum              INTEGER,
    ser_og_his_ne2_dist_A         REAL,
    his_nd1_asp_od_dist_A         REAL,
    ser_his_asp_angle_deg         REAL,
    oxyanion_n1_dist_A            REAL,
    oxyanion_n2_dist_A            REAL,
    cleft_width_A                 REAL,
    cleft_depth_A                 REAL,
    aromatic_clamp_residues_json  TEXT,
    solvent_accessibility_json    TEXT,
    FOREIGN KEY(candidate_id) REFERENCES structures(candidate_id)
);

-- ============================================================
-- Reference set: characterised PETase/MHETase-like plus ESTHER negatives
-- ============================================================
CREATE TABLE IF NOT EXISTS characterised_enzymes (
    enzyme_id                 TEXT PRIMARY KEY,   -- 'IsPETase', 'FAST-PETase', ...
    uniprot                   TEXT,
    pdb_ids_json              TEXT,
    organism                  TEXT,
    family                    TEXT,      -- petase_like | mhetase_like | cutinase | lipase | carboxylesterase | other
    sequence                  TEXT,
    seq_length                INTEGER,
    is_positive               INTEGER,   -- PAZy positive (0/1)
    is_negative               INTEGER,   -- ESTHER hard negative (0/1)
    is_near_miss              INTEGER,   -- active on soluble esters, not PET (0/1)
    esther_family             TEXT,
    matched_positive_id       TEXT,      -- for negatives: which positive it was matched to
    topt_c                    REAL,
    ph_opt                    REAL,
    activity_substrate_notes  TEXT,
    source_ref                TEXT,      -- PAZy | ESTHER | DOI
    pdb_release_date          TEXT,      -- for prospective-holdout cutoff logic
    added_at                  TEXT,
    FOREIGN KEY(matched_positive_id) REFERENCES characterised_enzymes(enzyme_id)
);

CREATE INDEX IF NOT EXISTS idx_char_family ON characterised_enzymes(family);
CREATE INDEX IF NOT EXISTS idx_char_pos    ON characterised_enzymes(is_positive, is_negative);

-- ============================================================
-- Harmonised literature activity data (spec section 5.4)
-- ============================================================
CREATE TABLE IF NOT EXISTS activity_measurements (
    id                     INTEGER PRIMARY KEY AUTOINCREMENT,
    enzyme_id              TEXT NOT NULL,
    substrate_form         TEXT,      -- amorphous_film | crystalline_powder | nanoparticle
    crystallinity_pct      REAL,
    temperature_c          REAL,
    ph                     REAL,
    buffer                 TEXT,
    assay_duration_h       REAL,
    product_measured       TEXT,      -- TPA_HPLC | turbidity | weight_loss
    rate_value             REAL,
    rate_units             TEXT,
    comparable_group_id    TEXT,      -- rows sharing this id ARE mutually comparable
    ordinal_rank_in_paper  INTEGER,   -- fallback when harmonisation fails (spec section 5.4)
    source_doi             TEXT,
    extracted_at           TEXT,
    extraction_confidence  TEXT,      -- high | medium | low
    FOREIGN KEY(enzyme_id) REFERENCES characterised_enzymes(enzyme_id)
);

-- ============================================================
-- Evaluation (spec section 8)
-- ============================================================
CREATE TABLE IF NOT EXISTS eval_splits (
    sequence_id            TEXT NOT NULL,   -- candidate_id or enzyme_id
    cluster_id             TEXT,
    identity_threshold     INTEGER,         -- 30 or 50
    split                  TEXT,            -- train | val | test
    is_prospective_holdout INTEGER,
    PRIMARY KEY(sequence_id, identity_threshold)
);

CREATE TABLE IF NOT EXISTS training_runs (
    run_id                     INTEGER PRIMARY KEY AUTOINCREMENT,
    head_name                  TEXT,     -- pet_activity (only head in v1)
    model_version              TEXT,
    n_positives                INTEGER,
    n_negatives                INTEGER,
    pu_class_prior             REAL,
    pu_method                  TEXT,     -- elkan_noto | bagging
    calibration_method         TEXT,     -- platt | isotonic
    cluster_identity_threshold INTEGER,
    auc                        REAL,
    average_precision          REAL,
    brier_score                REAL,
    retrieval_baseline_auc     REAL,     -- HMMER E-value rank: the bar to beat
    trained_at                 TEXT,
    config_json                TEXT
);

-- ============================================================
-- App-level (footer version display, lightweight visitor counter)
-- ============================================================
CREATE TABLE IF NOT EXISTS app_state (
    key    TEXT PRIMARY KEY,   -- 'current_data_version', 'current_model_version'
    value  TEXT
);

CREATE TABLE IF NOT EXISTS visits (
    ip_hash     TEXT PRIMARY KEY,
    first_seen  TEXT,
    last_seen   TEXT,
    hits        INTEGER DEFAULT 0
);
```

Deferred to v2 (schema sketch only, not built in v1): `variant_scores` (Engineer tab, position by amino-acid heatmap cache), `compare_sessions` (Compare tab), `submit_requests` (rate-limiting and audit for live `/submit`). These extend cleanly off `candidates` and `scores` without a v1 migration.

---

## 4. Offline compute time and disk footprint, with the arithmetic

**Metagenome collection (decision 1, sized against the real 62 GB free budget):** one MGnify plastisphere-or-landfill assembly analysis study. Assume, conservatively, ~3 million predicted protein sequences at ~300 aa average, ~400 bytes per sequence including FASTA header, giving **~1.2 GB** as protein FASTA. Confirm via the MGnify API's own size metadata *before* the bulk download starts (risk 5): breadth-first across MGnify plus JGI IMG/M plus OceanDNA plus Tara Oceans would be tens to hundreds of GB and simply does not fit.

**MMseqs2 prefilter:** ~200 seed sequences (positives) against ~3M targets at high sensitivity (`-s 7.5`). Published throughput at this scale is order 10^5 to 10^6 comparisons/sec/core; on the M1 Max's 10 cores, a defensible estimate is **2 to 4 hours**. Replace this with a 1,000-sequence smoke test in week 1 rather than trusting it outright.

**HMMER sensitive pass:** MMseqs2 typically narrows 3M to order 10,000 to 50,000 survivors. `hmmsearch` against ~15 profile HMMs (PETase family plus ESTHER families), at throughput on the order of tens of sequences/sec/core, scales to roughly **1 to 3 hours** on 10 cores.

**ESM-2 t12-35M embedding:** after the triad/oxyanion-hole filter, expect candidates in the low thousands (say 2,000 to 5,000) plus a few hundred positives and negatives, so ~5,000 to 10,000 sequences total. At a conservative CPU throughput of 5 to 10 seq/s for a 35M-param encoder-only model: 10,000 / 7 = **25 to 40 minutes**. Generously pessimistic at 2 seq/s it is still only **~1.4 hours**.

**ESMFold, top 50:** no directly comparable in-house number exists (BoltzMaker's timings are for Boltz-2's diffusion process, architecturally heavier). Estimate 30 s to 5 min per structure for a ~300-residue monomer, giving **50 x (0.5 to 5 min) = 25 min to 4.2 hours**. **Smoke-test structure 1 before committing to structure 50** (risk 4). If it is bad, the fallback is Boltz-2 via BoltzMaker mono-fold, for which real numbers exist on this exact machine (T4 lysozyme, 164 residues, monomer, one ligand: **3 min 16 s** including preflight and affinity prediction), so 50 structures at a similar order is **~2.7 to 4 hours**, tolerable as a one-off offline batch.

**Disk footprint, v1, summed:**

| Item | Arithmetic | Size |
|---|---|---|
| Metagenome protein FASTA | ~3M seq x 400 B | ~1.2 GB |
| MMseqs2 DB indices | ~1.5x FASTA | ~1.8 GB |
| HMMER/MMseqs2 hit tables | small text | <100 MB |
| ESM-2 t12-35M weights | fp32, 35M params | ~150 MB (cached once) |
| Candidate embedding cache | 5,000 x 480-dim x 4 B | ~9.6 MB |
| Top-50 mmCIF (predicted, no waters) | 50 x ~200 KB (AlphaFraud's experimental-PDB average is ~350 KB/file over ~2,985 files = 1.0 GB; a predicted monomer without waters or alt-confs runs lighter) | ~10 MB |
| PDB reference structures (~25 characterised enzymes) | 25 x ~400 KB | ~10 MB |
| SQLite catalogue (~5,000 rows including sequence text) | 5,000 x ~1.5 KB | ~7.5 MB |
| **Total** | | **~3.2 GB, comfortably inside the 62 GB local budget** |

The droplet side is trivial by comparison: the deployed SQLite file, top-50 mmCIF and a thin web venv (no torch) sum to well under 1 GB. No disk risk for v1.

---

## 5. Memory check: on-demand inference alongside gunicorn

**This could not be measured directly** (no SSH access to the droplet from the planning sandbox). The numbers below are principled estimates from known worker counts and library footprints, explicitly not measurements, and **task 8.4 makes the real `free -m` check the first action of deploy**, so this section gets replaced by ground truth before any v2 code is written.

Estimated current baseline resident memory across the five existing apps plus OS and nginx:

| App | Workers | Estimated RSS/worker | Subtotal |
|---|---|---|---|
| OS, sshd, systemd, nginx | | | ~200 MB |
| AlphaFraud (numpy/scipy/pandas/biotite/tmtools) | 3 sync | ~190 MB | ~570 MB |
| chem_sage-web (thin gevent, proxies to remote HF Space) | 1 gevent | ~220 MB | ~220 MB |
| chatPDB-web (same pattern) | 1 gevent | ~150 MB | ~150 MB |
| BoltzMaker-web (thin flask-only venv, subprocess for heavy work) | 3 sync | ~75 MB | ~225 MB |
| FlexAppeal (baseline, MemoryMax cap 1200 MB) | 2 sync | ~120 MB | ~240 MB |
| **Baseline total** | | | **~1.6 GB** |

Nominal headroom on 3.8 GB: **~2.2 GB**. But FlexAppeal can spike to its full 1200 MB cap and BoltzMaker's subprocess can spike to its 1 GB cap. A simultaneous worst case (rare, but the reason FlexAppeal's cgroup exists at all) could push resident use to **~3.3 GB**, leaving as little as **~500 MB** of true safety margin.

**v1's own footprint (no torch, no ESM-2 in the request path: v1 ships Home, Catalogue, Candidate and Methods only, all reading precomputed SQLite):** 1 to 2 sync workers x ~150 MB = **150 to 300 MB added**. Fits easily. This is why v1's memory story is genuinely not a problem.

**v2's footprint (Submit and Engineer, live ESM-2 CPU inference) is the real constraint, and the spec's "under 2 GB" undersells it:**

- torch CPU import alone: ~300 to 500 MB RSS (MKL/OpenMP shared libraries)
- transformers overhead: ~50 to 100 MB
- ESM-2 t12-35M weights resident: ~150 to 200 MB
- Per-request activation buffers (~500-residue sequence): ~50 to 150 MB
- **Total per worker with the model loaded: ~700 MB to 1 GB**

One such worker alone would consume most of the realistic 500 MB to 2 GB headroom above, before any concurrent load from sibling apps. Two such workers (a naive `workers=2`) would not fit at all.

**Recommendation for v2, documented now so v1's architecture does not foreclose it:**

1. `workers=1` for the PANTS gunicorn service, always. This is a low-traffic personal tool; serialising `/submit` requests is an acceptable trade for not duplicating a 700 MB to 1 GB process.
2. Lazy-import torch and transformers **inside** the `/submit` and `/engineer` view functions, not at module load, so v1's browsing workers never pay the cost and the cost is paid once per worker on first live request.
3. `MemoryMax` of roughly 900M to 1G on the PANTS web service, matching FlexAppeal's and BoltzMaker-web's precedent, purely as a backstop so a runaway PANTS request cannot OOM the other five apps.
4. **Fallback if the real numbers from 8.4 do not work:** offload ESM-2 inference to a remote HF Space, exactly as chem_sage-web and chatPDB-web already do (`chat_remote.py`'s thin-proxy pattern, near-zero droplet marginal memory). One caveat from chatPDB's own provisioning notes: ZeroGPU's anonymous quota is roughly 85 s/day and needs a Bearer token, and a non-ZeroGPU always-on small CPU Space would avoid the 60 to 120 s cold start chatPDB documents, which the spec's under-10-second budget could not otherwise absorb.

---

## 6. Section 10 decisions: which block v1

| # | Decision | v1-blocking? | Recommendation |
|---|---|---|---|
| 1 | Metagenome collection | **Yes** (v1 needs one collection) | Plastisphere/landfill, narrow, not breadth-first. The 62 GB real local budget rules out breadth-first outright; plastic-surface-associated communities have higher expected true-positive density, which matters for the free-baseline test in spec section 2; and it matches the narrow-then-expand pattern AlphaFraud itself used. |
| 2 | MHETase shared or parallel pipeline | No | Defer entirely to v2, as a parallel pipeline with its own seed and negatives. MHETase's tannase/feruloyl-esterase fold genuinely differs (spec section 2, point 5), so sharing the PETase HMM profile would miss it regardless of when it is built. |
| 3 | PU class prior | **Yes**: v1's single head cannot train without it | Elkan-Noto, not bagging PU: one hyperparameter, defensible at this sample size. Estimate a conservative literature-informed prior (order 1 to 5% of triad-complete candidates in a plastic-associated environment genuinely PET-active) and **sensitivity-test rankings at 1, 3, 5 and 10%** explicitly in the Methods tab, rather than presenting one unexamined number. |
| 4 | Regression head versus ranking-only | **Already resolved by spec section 13** | v1's scope line says "one calibrated activity head" (binary, per section 5.5's first row). Not open, decided. Stated here so it does not get re-litigated. |
| 5 | Structure budget (count, method, disk) | **Partially**: count is fixed at 50 by section 13, method needs a v1 decision | ESMFold only for v1: no ligand co-fold is needed (Engineer, which would want it, is v2), and ESMFold is a single forward pass versus Boltz-2's diffusion, which is minutes-to-tens-of-minutes heavier per structure on this hardware per BoltzMaker's real timings. Reserve Boltz-2/BoltzMaker for v2's ligand overlay. |
| 6 | Submit/Engineer public or gated | No | v1 ships neither route. Decide once section 5's real memory numbers exist. Leaning public-with-cgroup-containment if the numbers work, HF-Space offload if they do not. |
| 7 | Licensing and dataset release | No | Policy question, not a v1 build blocker. Populate `data_sources.license` from v1 onward (already in the schema) so the decision does not require re-deriving provenance later. |

---

## 7. v1 scope line

**Ships:**

- Repo scaffold (spec section 11 layout), manifest infrastructure, `PANTS.py` CLI (`init`/`run`/`serve`, mirroring `AlphaFraud.py`)
- One metagenome collection (plastisphere/landfill), MMseqs2 prefilter plus HMMER sensitive pass, triad/oxyanion-hole completeness filter, discard counts reported
- ESTHER hard negative set (length, identity and taxonomy matched, near-misses included), sanity-checked against a trivial baseline before any embedding work (risk 1)
- Characterised reference set: FAST-PETase, IsPETase, DuraPETase, HotPETase, LCC-ICCG, Z1-PETase, TurboPETase, Cut190\*\*SS, HGMP01, MHETase, plus PAZy's own list, explicitly not claimed exhaustive
- ESM-2 t12-35M frozen embeddings (offline, CPU)
- **One** calibrated PET activity head: PU-corrected (Elkan-Noto), Platt/isotonic calibration, class-prior sensitivity analysis
- Full evaluation protocol from spec section 8 (cluster splits, leave-one-family-out, retrieval baseline, calibration diagram, prospective holdout)
- Structures for the top 50 candidates: ESMFold, superposed onto IsPETase 6EQE (TM-score and RMSD), active-site geometry (triad, oxyanion hole, cleft width and depth, aromatic clamp, solvent accessibility)
- Home, Catalogue, Candidate and Methods tabs, live at `pants.mdeller.com`, deploy mirroring AlphaFraud exactly (gunicorn, nginx, systemd, certbot, port 8005), added to `mdeller-landing/apps.json`
- Full SQLite schema including manifest and versioning tables

**Home tab hero plot, resolved honestly (not spec-literal, since there is no Topt head):** activity-probability distribution across the catalogue with characterised enzymes overlaid as anchors, or activity probability versus identity-to-nearest-characterised-enzyme. The latter directly visualises the spec section 2 thesis: detection is solved, discrimination is not. Label it explicitly in the UI as the v1 placeholder for the spec'd activity-versus-Topt plot, which ships once the Topt head lands.

**Deferred to v2:** MHET activity head and its parallel pipeline; Topt, expressibility and aggregation heads (so the therapeutic composite becomes meaningfully multi-axis only in v2); Boltz-2/BoltzMaker ligand co-folds; structure budget beyond 50; Compare, Engineer, Submit and Analysis tabs entirely; additional metagenome collections; the dataset-licensing decision.

---

## 8. Risks, ranked, with mitigations

1. **Hard negative construction is the single point of failure (the spec's own words).** A poorly matched negative set lets the model win on a shortcut (length, taxonomic origin) rather than the real discriminative signal. *Mitigation:* before touching ESM-2, run a trivial baseline (logistic regression on amino-acid composition or raw length alone) on the matched set. If that already separates positives from negatives well, the matching failed and needs redoing. One afternoon, before Phase 4.

2. **Cluster-split leakage would make "calibrated" a lie.** With order 10^2 positives, a wrongly clustered split (or a near-miss cutinase leaking across train and test) would produce evaluation numbers that look fine and are not. *Mitigation:* one assertion, that max pairwise identity between train and test is at or below the declared threshold, run at Phase 5 kickoff. Minutes, with mmseqs2 already installed.

3. **PAZy curation has no API, and "and ALL others you can find" is open-ended.** Manual curation is slower than day-count estimates usually assume. *Mitigation:* timebox to the named list plus one pass of PAZy's own table; treat "ALL others" as an ever-expanding Methods-tab table, not a v1 blocker. Spend 2 hours on day 1, measure rows per hour, extrapolate before trusting Phase 1's schedule.

4. **ESMFold on Apple Silicon is asserted, not verified.** HF's reimplementation should avoid the CUDA-only `openfold` build, but BoltzMaker's README documents real, specific MPS pain (a hard crash above ~1250 residues, a silent CPU fallback for `torch.linalg.svd`) for a different but related model. *Mitigation:* one structure, timed, checked with gemmi or biotite parse, before committing to all 50. Fallback with real numbers already in hand: Boltz-2 via BoltzMaker mono-fold, ~2.7 to 4 hours for 50 structures.

5. **62 GB free local disk, not "plenty".** Downloading the wrong-sized metagenome collection could exhaust it mid-run, and a misdirected download into `~/Documents` risks iCloud eviction and thrashing. *Mitigation:* set up the symlink pattern before any download (Phase 0.4); query the chosen MGnify study's size via its own API metadata before the bulk pull, and refuse anything over ~30 GB. Minutes, as the first action of Phase 1.4.

6. **The droplet memory numbers in section 5 are estimates, not measurements.** If real headroom is smaller than estimated, v1 (no torch) is still safe, but v2's rollout risk transfers entirely to whichever app is unlucky enough to be under load at the same moment. *Mitigation:* `free -m` and `systemctl status` on all five live services, as literally the first task of Phase 8, before writing any v2 memory-sensitive code.

7. **Active-site geometry code is genuinely new, with nothing to port.** Getting cleft width or the aromatic clamp wrong would quietly undermine exactly the feature spec section 6 flags as most discriminative. *Mitigation:* validate against known PDB structures (IsPETase versus a confirmed non-PET-active cutinase or lipase) before running on any predicted structure. Half a day, doable in parallel with Phases 2 to 5.

8. **Python 3.14 and torch compatibility is unconfirmed against this house's own precedent.** chem_sage deliberately used 3.11+ for its torch stack and the reason is not documented. *Mitigation:* `pip install torch` in a scratch 3.14 venv as literally the first command of Phase 0; fall back to 3.11/3.12 (chem_sage's proven choice) if it fails.

---

## 9. Reference files to port from

PANTS itself has only `pants_spec_v1.md` today, so every pattern below comes from a sibling project already running in production.

- `AlphaFraud/alphafraud/db.py`: SQLite connect/WAL/manifest pattern to port
- `AlphaFraud/alphafraud/compare.py`: TM-align/tmtools superposition, including the `_TMALIGN_LOCK` thread-safety guard, to reuse for task 6.3
- `AlphaFraud/alphafraud/http.py`: requests + tenacity retry wrapper for tasks 1.2 and 1.3
- `AlphaFraud/alphafraud/pdb.py`: RCSB GraphQL fetch for task 1.5
- `AlphaFraud/deploy/`: gunicorn.conf.py, nginx conf, systemd units, provision.sh, deploy.sh, the exact deploy skeleton to fork onto port 8005
- `AlphaFraud/alphafraud/static/brand.css`: house theme for task 7.1
- `FlexAppeal/deploy/flexappeal-web.service`: MemoryMax cgroup pattern for v2 inference containment
- `mdeller-landing/apps.json`: the exact file to add the `pants` entry to
- `TopPDBLX/PHASE0_PLAN.md`: iCloud-avoidance symlink pattern and manifest convention
