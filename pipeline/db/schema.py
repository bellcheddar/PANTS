"""The PANTS SQLite schema, as one executescript-able string.

Conventions ported from AlphaFraud's db.py, which has been in production since 2026-07-14:

  - Large JSON blobs (per-residue arrays, position lists) live in sidecar tables or JSON
    columns, never in the hot path. The handful of columns used for sorting, filtering and
    the Catalogue table are promoted to real columns with indices on them. AlphaFraud
    learned this the hard way: an uncovered index on a 28-column table cost 62 ms per page
    load until a covering index fixed it.
  - Every pipeline stage opens a row in `runs` and writes a row in `manifests` recording
    input/output hashes, tool versions and the git commit (spec section 12).
  - v2 columns (mhet_activity_prob, topt_pred_c, expressibility_prob, aggregation_score)
    exist now and stay NULL through v1, so adding the v2 heads needs no migration.
"""

from __future__ import annotations

# Bump when SCHEMA changes in a way that invalidates an existing pants.db. Recorded in
# every manifest so a stale database is obvious rather than silently mis-read.
SCHEMA_VERSION = 16

SCHEMA = """
-- ============================================================
-- Provenance / manifest / versioning (spec section 12)
-- ============================================================
CREATE TABLE IF NOT EXISTS runs (
    id           INTEGER PRIMARY KEY AUTOINCREMENT,
    stage        TEXT NOT NULL,          -- recall | negatives | embed | train | structure | harmonise
    label        TEXT NOT NULL,          -- run date or dataset version tag
    started_at   TEXT,
    finished_at  TEXT,
    n_input      INTEGER DEFAULT 0,
    n_output     INTEGER DEFAULT 0,
    n_discarded  INTEGER DEFAULT 0,      -- e.g. candidates failing the triad filter
    status       TEXT DEFAULT 'running', -- running | done | error
    params_json  TEXT
);

CREATE TABLE IF NOT EXISTS manifests (
    id                  INTEGER PRIMARY KEY AUTOINCREMENT,
    run_id              INTEGER,
    stage               TEXT NOT NULL,
    input_hashes_json   TEXT,            -- {path: sha256}
    output_hashes_json  TEXT,
    tool_versions_json  TEXT,            -- {hmmer: '3.4', mmseqs2: '18-8cc5c', ...}
    model_version       TEXT,
    schema_version      INTEGER,
    git_commit          TEXT,
    wall_time_s         REAL,
    written_at          TEXT,
    FOREIGN KEY(run_id) REFERENCES runs(id)
);

CREATE TABLE IF NOT EXISTS data_sources (
    name          TEXT PRIMARY KEY,      -- PAZy | ESTHER | PDB | MGnify:<study_id> | Meltome
    version       TEXT,
    retrieved_at  TEXT,
    n_records     INTEGER,
    license       TEXT,                  -- populated from v1 so decision 7 needs no re-derivation
    source_url    TEXT
);

-- ============================================================
-- Candidates (recall output, catalogue backbone)
-- ============================================================
CREATE TABLE IF NOT EXISTS candidates (
    candidate_id                   TEXT PRIMARY KEY,  -- 'PANTS-' + first 12 hex of SHA256(sequence)
    source_environment             TEXT,              -- plastisphere | landfill | marine | wastewater | compost
    assembly_id                    TEXT,
    contig_id                      TEXT,
    sample_accession               TEXT,
    taxon_name                     TEXT,
    ncbi_taxid                     INTEGER,
    taxonomy_lineage               TEXT,
    sequence                       TEXT NOT NULL,
    seq_length                     INTEGER,
    has_complete_triad             INTEGER,           -- 0/1 filter gate, spec section 5.1
    triad_positions_json           TEXT,
    oxyanion_hole_positions_json   TEXT,
    recall_method                  TEXT,              -- mmseqs2 | hmmer
    recall_evalue                  REAL,
    recall_bitscore                REAL,
    recall_profile_identity        REAL,
    nearest_characterised_id       TEXT,
    nearest_characterised_identity REAL,
    discovered_run_id              INTEGER,
    first_seen_at                  TEXT,
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

-- Cleft width and the aromatic clamp are the features most likely to separate real
-- polyesterases from soluble-ester-only esterases (spec section 6), so they are promoted
-- to real columns and surfaced prominently rather than buried in JSON.
CREATE TABLE IF NOT EXISTS geometry (
    candidate_id                  TEXT PRIMARY KEY,
    triad_ser_resnum              INTEGER,
    triad_his_resnum              INTEGER,
    triad_asp_resnum              INTEGER,
    ser_og_his_ne2_dist_A         REAL,
    his_nd1_asp_od_dist_A         REAL,
    ser_his_asp_angle_deg         REAL,
    oxyanion_n1_dist_A            REAL,
    oxyanion_n1_resnum            INTEGER,
    oxyanion_n2_resnum            INTEGER,
    oxyanion_n2_angle_deg         REAL,
    oxyanion_n2_dist_A            REAL,
    cleft_width_A                 REAL,
    cleft_depth_A                 REAL,
    aromatic_clamp_residues_json  TEXT,
    solvent_accessibility_json    TEXT,
    FOREIGN KEY(candidate_id) REFERENCES structures(candidate_id)
);

CREATE INDEX IF NOT EXISTS idx_geometry_cleft ON geometry(cleft_width_A);

-- Structures and geometry for the CHARACTERISED enzymes, mirroring the candidate tables.
--
-- Separate tables rather than a shared one because `structures` and `geometry` are keyed
-- on candidate_id with a foreign key into `candidates`, and a characterised enzyme is not
-- a candidate. Widening those keys would let a metagenomic hit and a reference enzyme
-- collide in the same namespace, which is precisely the distinction the whole evaluation
-- rests on.
--
-- `source` records where the coordinates came from, because it is not one thing: an
-- experimental PDB entry where one exists, otherwise ESMFold for the engineered variants,
-- which have neither a PDB entry nor a UniProt accession of their own. A geometric
-- comparison across a mixture of crystal structures and predictions is only honest if the
-- reader can see which is which.
CREATE TABLE IF NOT EXISTS reference_structures (
    enzyme_id               TEXT PRIMARY KEY,
    source                  TEXT,        -- 'pdb' | 'alphafold' | 'esmfold'
    source_id               TEXT,        -- PDB id or UniProt accession
    coord_path              TEXT,        -- viewer PDB, superposed, under app/static/reference_structures
    plddt_mean              REAL,        -- predictions only; NULL for experimental
    resolution_A            REAL,        -- experimental only
    rmsd_ca_to_ispetase_A   REAL,
    superposition_reference TEXT DEFAULT '6EQE',
    n_residues              INTEGER,
    built_at                TEXT,
    model_version           TEXT,
    FOREIGN KEY(enzyme_id) REFERENCES characterised_enzymes(enzyme_id)
);

CREATE TABLE IF NOT EXISTS reference_geometry (
    enzyme_id                     TEXT PRIMARY KEY,
    triad_ser_resnum              INTEGER,
    triad_his_resnum              INTEGER,
    triad_asp_resnum              INTEGER,
    ser_og_his_ne2_dist_A         REAL,
    his_nd1_asp_od_dist_A         REAL,
    ser_his_asp_angle_deg         REAL,
    oxyanion_n1_dist_A            REAL,
    oxyanion_n1_resnum            INTEGER,
    oxyanion_n2_dist_A            REAL,
    oxyanion_n2_resnum            INTEGER,
    oxyanion_n2_angle_deg         REAL,
    cleft_width_A                 REAL,
    cleft_depth_A                 REAL,
    n_cleft_residues              INTEGER,
    aromatic_clamp_residues_json  TEXT,
    FOREIGN KEY(enzyme_id) REFERENCES reference_structures(enzyme_id)
);

-- ============================================================
-- Reference set: characterised PETase/MHETase-like plus ESTHER hard negatives
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
    is_near_miss              INTEGER,   -- active on soluble esters, NOT on PET (0/1).
                                         -- These define the decision boundary (spec 5.2).
    esther_family             TEXT,
    taxonomy_lineage          TEXT,      -- full lineage; phylum is used to match negatives
                                         -- to positives on composition-driving taxonomy
    matched_positive_id       TEXT,      -- for negatives: which positive it was matched to
    topt_c                    REAL,
    ph_opt                    REAL,
    activity_substrate_notes  TEXT,
    source_ref                TEXT,      -- PAZy | ESTHER | DOI
    pdb_release_date          TEXT,      -- prospective-holdout cutoff logic (spec section 8)
    is_fragment               INTEGER,   -- UniProt's own Fragment flag, not a length guess
    excluded_from_training    INTEGER DEFAULT 0,  -- marked, never deleted: the catalogue
    exclusion_reason          TEXT,               -- stays complete and the exclusion is
                                                  -- auditable rather than invisible
    added_at                  TEXT,
    FOREIGN KEY(matched_positive_id) REFERENCES characterised_enzymes(enzyme_id)
);

CREATE INDEX IF NOT EXISTS idx_char_family ON characterised_enzymes(family);
CREATE INDEX IF NOT EXISTS idx_char_pos    ON characterised_enzymes(is_positive, is_negative);

-- ============================================================
-- Harmonised literature activity data (spec section 5.4)
-- ============================================================
-- Published activity numbers are NOT comparable across papers. Rows sharing a
-- comparable_group_id are mutually comparable; where harmonisation is impossible,
-- ordinal_rank_in_paper carries the within-paper ranking instead and rate_value stays NULL.
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
    parameter_type         TEXT,      -- km | kcat | vmax | topt | ph_opt | specific_activity
    rate_value             REAL,
    rate_units             TEXT,
    raw_text               TEXT,      -- the source statement, kept verbatim. Optima are
                                      -- reported as prose ("Optimum pH is 8.5 with
                                      -- pNP-butyrate"), and the parsed number alone loses
                                      -- the substrate and the conditions it depends on.
    evidence_code          TEXT,      -- ECO:0000269 = experimental evidence from a
                                      -- publication. Anything weaker is not a measurement.
    comparable_group_id    TEXT,
    ordinal_rank_in_paper  INTEGER,
    source_doi             TEXT,
    extracted_at           TEXT,
    extraction_confidence  TEXT,      -- high | medium | low
    FOREIGN KEY(enzyme_id) REFERENCES characterised_enzymes(enzyme_id)
);

CREATE INDEX IF NOT EXISTS idx_activity_enzyme ON activity_measurements(enzyme_id);

-- ============================================================
-- Evaluation (spec section 8)
-- ============================================================
-- Splits are by sequence CLUSTER at 30% and 50% identity, never by sequence. Both
-- thresholds are stored side by side, hence the composite primary key.
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
    head_name                  TEXT,     -- pet_activity (the only head in v1)
    model_version              TEXT,
    n_positives                INTEGER,
    n_negatives                INTEGER,
    pu_class_prior             REAL,
    pu_method                  TEXT,     -- elkan_noto | bagging
    calibration_method         TEXT,     -- platt | isotonic
    cluster_identity_threshold INTEGER,
    auc                        REAL,
    average_precision          REAL,
    brier_score                REAL,     -- calibration is the point, so this is not optional
    retrieval_baseline_auc     REAL,     -- HMMER E-value rank: the bar the model must beat
    composition_baseline_auc   REAL,     -- amino-acid composition + length, cluster-grouped.
                                         -- Reported ALONGSIDE the model score permanently,
                                         -- not just used as a pre-training gate: it sat at
                                         -- 0.84 after five matching axes, so any model
                                         -- claim has to clear it as well as the E-value
                                         -- baseline (PHASE1_FINDINGS.md).
    n_positive_clusters        INTEGER,  -- independent units, NOT the raw positive count
    evidence_level             TEXT,     -- protein-evidence | predicted | mixed
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
"""


# Additive column migrations, applied by init_schema().
#
# The schema is CREATE TABLE IF NOT EXISTS throughout, which means a new column in the
# definition above NEVER reaches a database that already exists: the CREATE is skipped and
# the column silently stays missing. Every deployment of this project so far has been one
# such database, so a column added without an entry here is a column that exists only on a
# fresh machine.
#
# Each entry is (table, column, declaration) and is applied only when the column is
# absent, so this is safe to run on every startup.
COLUMN_MIGRATIONS = [
    # v6: oxyanion hole detection now records WHICH residues it identified, not just
    # distances, because the previous version reported a distance to the wrong residue and
    # nothing in the stored output could reveal that.
    ("geometry", "oxyanion_n1_resnum", "INTEGER"),
    ("geometry", "oxyanion_n2_resnum", "INTEGER"),
    ("geometry", "oxyanion_n2_angle_deg", "REAL"),
    # Which family definition admitted a within-family negative: cluster | profile | both.
    ("characterised_enzymes", "within_family_basis", "TEXT"),
    # v8: identity to the wild type at the root of this enzyme's lineage, and which enzyme
    # that is. Precomputed because the alignment needs biotite and the web virtual
    # environment is deliberately thin.
    ("characterised_enzymes", "lineage_wt_id", "TEXT"),
    ("characterised_enzymes", "identity_to_lineage_wt", "REAL"),
    # v9: one-line description for the summary table. A curated condensation of the notes
    # this project already holds, kept OUT of activity_measurements because most of these
    # describe what an enzyme IS rather than reporting a measured value.
    ("characterised_enzymes", "headline", "TEXT"),
    # v13: the name the literature uses. The bulk sets are keyed by accession -- PAZy:1,
    # ESTHER:C9ZCR8 -- which is stable and unreadable. PAZy:1 is IsPETase.
    ("characterised_enzymes", "common_name", "TEXT"),
    # v10: how the structure file's residue numbering relates to the stored sequence.
    # structure_resnum = sequence_position + seq_offset. A crystallised construct is
    # numbered by its depositors and need not agree with a precursor sequence: three of
    # twenty references were already out by +42, +26 and -40, which silently made the
    # sequence panel mark the wrong residues as catalytic.
    ("reference_structures", "seq_offset", "INTEGER"),
    # v11: for each substituted residue, how far it sits from the catalytic machinery.
    # Precomputed from the superposed coordinates: measuring it at request time would put
    # geometry in the web process, which is deliberately thin.
    ("reference_structures", "mutation_geometry_json", "TEXT"),
    # v12: catalytic residues the deposit has lost. Non-empty means the coordinates are an
    # INACTIVATED crystallisation construct rather than the working enzyme.
    ("reference_structures", "catalytic_knockout_json", "TEXT"),
    # v14: the date UniProt first made the entry public. The prospective holdout in the
    # evaluation protocol needs a "what was known when" axis, and `pdb_release_date` is
    # empty for the whole catalogue because most of these enzymes have no deposit. The
    # UniProt date is the earliest defensible proxy for when the sequence entered the
    # public record, and it exists for anything with an accession.
    ("characterised_enzymes", "uniprot_first_public", "TEXT"),
    # v15: the paper an enzyme's activity was reported in, as a column rather than as
    # free text. The PAZy import wrote "Primary reference doi:..." into
    # activity_substrate_notes, which preserved the citation but made it unqueryable --
    # and the whole ordinal-ranking route depends on GROUPING enzymes by paper, since a
    # ranking only exists where one protocol assayed several enzymes.
    ("characterised_enzymes", "primary_doi", "TEXT"),
    # v16: the sequence cluster an enzyme belongs to in its source paper's own framework,
    # which is what the paper's ecological metadata is keyed on.
    ("characterised_enzymes", "science_cluster", "INTEGER"),
]

# Tables added after the original schema. Same reason as COLUMN_MIGRATIONS: the CREATE
# statements above are skipped wholesale on an existing database.
LATE_TABLES = """
-- Ecological context for a cluster of related sequences, as reported by the paper that
-- defined the cluster. Kept at CLUSTER level because that is the level it was reported at:
-- attaching an isolation temperature to an individual enzyme would invent a precision the
-- source does not have.
CREATE TABLE IF NOT EXISTS sequence_clusters (
    source              TEXT NOT NULL,
    cluster_id          TEXT NOT NULL,
    n_members           INTEGER,
    isolation_sources_json TEXT,
    biomes_json         TEXT,
    habitats_json       TEXT,
    locations_json      TEXT,
    temperatures_json   TEXT,   -- raw strings, kept verbatim: "4,0 degree of C", "65 degrees F"
    temperature_median_c REAL,  -- parsed where parseable, NULL where not
    added_at            TEXT,
    PRIMARY KEY (source, cluster_id)
);

-- Sequences with no activity label at all. Deliberately NOT in characterised_enzymes:
-- nothing here has been characterised, and letting 23,000 unlabelled homologues into a
-- table whose counts are quoted as evidence would corrupt every total on the site.
CREATE TABLE IF NOT EXISTS unlabelled_sequences (
    seq_id        TEXT PRIMARY KEY,
    accession     TEXT,
    description   TEXT,
    organism      TEXT,
    sequence      TEXT NOT NULL,
    seq_length    INTEGER,
    source_ref    TEXT,
    source_doi    TEXT,
    added_at      TEXT
);
CREATE INDEX IF NOT EXISTS idx_unlabelled_len ON unlabelled_sequences(seq_length);
"""
