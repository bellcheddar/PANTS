"""Central configuration for PANTS.

All paths are anchored to the repository root (the parent of this package) so the code
runs identically from a laptop checkout and from /opt/pants on the droplet. Secrets and
the deploy target come from a gitignored .env (see .env.example); everything else has a
sensible default here.

This module is imported by BOTH the pipeline venv and the web venv, so it must never
import torch, transformers or anything else heavy. Keep it to the standard library plus
python-dotenv.
"""

from __future__ import annotations

import os
from pathlib import Path

try:
    from dotenv import load_dotenv
except ImportError:  # dotenv is optional at import time (e.g. bare python3 running `init`)
    def load_dotenv(*_args, **_kwargs):  # type: ignore
        return False

# --------------------------------------------------------------------------------------
# Paths
# --------------------------------------------------------------------------------------
PACKAGE_DIR = Path(__file__).resolve().parent
ROOT_DIR = PACKAGE_DIR.parent
DATA_DIR = ROOT_DIR / "data"

# data/raw and data/interim are SYMLINKS to ~/PANTSData, deliberately outside the
# iCloud-synced Documents tree: macOS "Optimize Mac Storage" evicts large files, and this
# machine has only ~62 GB free. Only data/processed (small, derived) lives in the repo.
RAW_DIR = DATA_DIR / "raw"
INTERIM_DIR = DATA_DIR / "interim"
PROCESSED_DIR = DATA_DIR / "processed"

MANIFEST_DIR = ROOT_DIR / "manifests"
MODELS_DIR = ROOT_DIR / "models"
EVAL_DIR = ROOT_DIR / "eval"
DB_PATH = Path(os.environ.get("PANTS_DB", ROOT_DIR / "pants.db"))

APP_DIR = ROOT_DIR / "app"
STATIC_DIR = APP_DIR / "static"
STRUCTURE_DIR = STATIC_DIR / "structures"   # precomputed mmCIF assets served statically

load_dotenv(ROOT_DIR / ".env")

# --------------------------------------------------------------------------------------
# Scientific constants
# --------------------------------------------------------------------------------------
# The therapeutic window PANTS ranks towards, versus the industrial optimum the published
# field targets (roughly 65 to 75 C, pH 8 to 9). See spec section 1.1.
THERAPEUTIC_TEMP_C = 37.0
THERAPEUTIC_PH = 7.4
THERAPEUTIC_TEMP_WINDOW_C = (30.0, 45.0)    # shaded band on the Home hero plot

# Superposition reference: IsPETase. 6EQE is the reference structure the geometry stage
# aligns every candidate onto (spec section 6.2).
ISPETASE_REFERENCE_PDB = "6EQE"

# Recall stage: candidates must carry a complete Ser-His-Asp triad to be scored at all.
# The count discarded here is reported (spec section 5.1).
REQUIRE_COMPLETE_TRIAD = True

# Sequence-cluster identity thresholds for evaluation splits. Both are reported, never
# a split by sequence (spec section 8).
CLUSTER_IDENTITY_THRESHOLDS = (30, 50)

# PU learning: the class prior is estimated, not assumed, and rankings are sensitivity
# tested across this grid with the result shown in the Methods tab (plan section 6).
PU_CLASS_PRIOR_DEFAULT = 0.03
PU_CLASS_PRIOR_GRID = (0.01, 0.03, 0.05, 0.10)

# Structure budget for v1: the top N candidates by activity probability get ESMFold.
STRUCTURE_BUDGET = int(os.environ.get("STRUCTURE_BUDGET", "50"))

# --------------------------------------------------------------------------------------
# Models
# --------------------------------------------------------------------------------------
# t12-35M, not 650M: ~150 MB fp32, the only ESM-2 that fits alongside gunicorn on a
# 3.8 GB droplet shared with five other apps (plan section 5).
ESM_MODEL = os.environ.get("ESM_MODEL", "facebook/esm2_t12_35M_UR50D")
ESM_EMBED_DIM = 480
ESMFOLD_MODEL = os.environ.get("ESMFOLD_MODEL", "facebook/esmfold_v1")

# --------------------------------------------------------------------------------------
# External tools (shelled out to, per the BoltzMaker house style)
# --------------------------------------------------------------------------------------
MMSEQS_BIN = os.environ.get("MMSEQS_BIN", "mmseqs")
HMMBUILD_BIN = os.environ.get("HMMBUILD_BIN", "hmmbuild")
HMMSEARCH_BIN = os.environ.get("HMMSEARCH_BIN", "hmmsearch")
MMSEQS_SENSITIVITY = float(os.environ.get("MMSEQS_SENSITIVITY", "7.5"))

# --------------------------------------------------------------------------------------
# External API endpoints
# --------------------------------------------------------------------------------------
UNIPROT_REST_URL = "https://rest.uniprot.org/uniprotkb/search"
RCSB_GRAPHQL_URL = "https://data.rcsb.org/graphql"
MGNIFY_API_URL = "https://www.ebi.ac.uk/metagenomics/api/v1"
HTTP_TIMEOUT = 60          # seconds per request
HTTP_MAX_RETRIES = 5

# Refuse a metagenome download above this size: local free disk is ~62 GB and the whole
# v1 footprint is budgeted at ~3.2 GB (plan risk 5).
MAX_COLLECTION_GB = float(os.environ.get("MAX_COLLECTION_GB", "30"))

# --------------------------------------------------------------------------------------
# Deploy / serving (from .env)
# --------------------------------------------------------------------------------------
DROPLET_SSH = os.environ.get("DROPLET_SSH", "")
DROPLET_PATH = os.environ.get("DROPLET_PATH", "/opt/pants")
SERVER_NAME = os.environ.get("SERVER_NAME", "pants.mdeller.com")
# Port 8005: AlphaFraud 8000, chem_sage-web 8001, chatPDB-web 8002, BoltzMaker-web 8003,
# FlexAppeal 8004. PANTS is app number six on the droplet.
BIND_ADDR = os.environ.get("BIND_ADDR", "127.0.0.1:8005")

# Data version shown in the app footer (spec section 12).
DATA_VERSION = os.environ.get("DATA_VERSION", "dev")


def ensure_dirs() -> None:
    """Create the runtime directory tree if missing. Safe to call repeatedly.

    RAW_DIR and INTERIM_DIR are symlinks to ~/PANTSData; mkdir on an existing symlink to a
    real directory is a no-op, and if the symlink is missing entirely we create a real
    directory rather than failing (the `init` command warns about that case).
    """
    for d in (DATA_DIR, PROCESSED_DIR, MANIFEST_DIR, MODELS_DIR, EVAL_DIR,
              STATIC_DIR, STRUCTURE_DIR):
        d.mkdir(parents=True, exist_ok=True)


def data_dirs_are_external() -> bool:
    """True when raw/ and interim/ really are symlinks out of the iCloud tree."""
    return RAW_DIR.is_symlink() and INTERIM_DIR.is_symlink()


# --------------------------------------------------------------------------------------
# Evidence tiers
#
# Which source_ref values denote a MEASURED positive, as opposed to one carrying automatic
# EC annotation or bare family membership. This is the axis the headline metric is
# reported on, so it must mean the same thing in the pipeline, the release bundle and the
# web app.
#
# It lives here because it did not, once: the app and the README disagreed (12 measured
# against 17) purely because two copies of this tuple had drifted apart. Copying it a
# third time would be the same bug waiting.
#
# The variant references are DERIVED from VARIANTS rather than typed out, so confirming a
# new mutation set cannot silently drop that variant out of the measured count. Each is a
# measured positive on the same criterion PAZy uses: someone built the enzyme, assayed it
# on PET and published the result.
# --------------------------------------------------------------------------------------
_BASE_MEASURED_TIERS = (
    "EC-experimental", "UniProt", "HGMP-measured", "PAZy-measured",
    # A construct taken from a crystal structure: it was expressed, crystallised and assayed.
    "PDB-construct",
)


def measured_tiers() -> tuple:
    """Every source_ref that counts as an experimentally measured positive."""
    from .recall.seeds import VARIANTS
    return _BASE_MEASURED_TIERS + tuple(sorted({v.reference for v in VARIANTS}))


MEASURED_TIERS = measured_tiers()
