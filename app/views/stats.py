"""Live instrumentation: what this project has processed, produced and is running on.

Every number here is read from the database or the host at request time. Nothing is cached
and nothing is hardcoded, so a figure that stops moving has stopped because the pipeline
stopped, not because a template went stale. That is the point of the page: it is the one
place where a claim about this project can be checked against the thing itself.

Stdlib only. The web virtual environment carries Flask and gunicorn and nothing else,
because the droplet has 3.9 GB shared with five other applications, so host metrics come
from `shutil.disk_usage`, `os.getloadavg` and `/proc` rather than from psutil.
"""

from __future__ import annotations

import json
import os
import platform
import shutil
import sqlite3
import threading
import time
from typing import Any, Dict, List, Optional

from flask import Blueprint, jsonify, render_template

from pipeline import config
from pipeline.db import connect
from pipeline.db.schema import SCHEMA_VERSION

bp = Blueprint("stats", __name__)

ENV_LABEL = {
    "compost": "Compost", "marine_plastisphere": "Marine plastisphere",
    "landfill": "Landfill", "wastewater": "Wastewater", "human_gut": "Human gut",
    "unknown": "Unknown", None: "Unspecified",
}

_STARTED = time.time()


def _one(conn, sql: str, params=()) -> Any:
    try:
        row = conn.execute(sql, params).fetchone()
        return row[0] if row else None
    except sqlite3.Error:
        return None


def _rows(conn, sql: str, params=()) -> List[Dict[str, Any]]:
    try:
        return [dict(r) for r in conn.execute(sql, params)]
    except sqlite3.Error:
        return []


def host_health() -> Dict[str, Any]:
    """Disk, memory and load, from the standard library.

    /proc exists on the droplet and not on macOS, so memory is reported where it can be
    read and omitted where it cannot, rather than guessed at from a different source that
    would not mean the same thing.
    """
    out: Dict[str, Any] = {}
    try:
        du = shutil.disk_usage(str(config.ROOT_DIR))
        out["disk"] = {"total_gb": round(du.total / 1e9, 1),
                       "used_gb": round(du.used / 1e9, 1),
                       "free_gb": round(du.free / 1e9, 1),
                       "used_pct": round(100 * du.used / du.total, 1)}
    except OSError:
        pass
    try:
        l1, l5, l15 = os.getloadavg()
        cpus = os.cpu_count() or 1
        out["load"] = {"1m": round(l1, 2), "5m": round(l5, 2), "15m": round(l15, 2),
                       "cpus": cpus, "pct_of_capacity": round(100 * l1 / cpus, 0)}
    except OSError:
        pass
    try:
        mem = {}
        with open("/proc/meminfo") as fh:
            for line in fh:
                k, _, v = line.partition(":")
                mem[k] = int(v.split()[0]) * 1024
        total, avail = mem.get("MemTotal", 0), mem.get("MemAvailable", 0)
        if total:
            out["memory"] = {"total_gb": round(total / 1e9, 2),
                             "available_gb": round(avail / 1e9, 2),
                             "used_pct": round(100 * (total - avail) / total, 1)}
        swap_t, swap_f = mem.get("SwapTotal", 0), mem.get("SwapFree", 0)
        if swap_t:
            out["swap"] = {"total_gb": round(swap_t / 1e9, 2),
                           "used_gb": round((swap_t - swap_f) / 1e9, 2)}
    except (OSError, ValueError):
        pass
    try:
        with open("/proc/uptime") as fh:
            out["host_uptime_h"] = round(float(fh.read().split()[0]) / 3600, 1)
    except (OSError, ValueError):
        pass
    out["process_uptime_h"] = round((time.time() - _STARTED) / 3600, 2)
    out["platform"] = f"{platform.system()} {platform.machine()}"
    out["python"] = platform.python_version()
    try:
        db = config.ROOT_DIR / "pants.db"
        out["database_mb"] = round(db.stat().st_size / 1e6, 1)
        out["database_age_min"] = round((time.time() - db.stat().st_mtime) / 60, 1)
    except OSError:
        pass
    out.update(_structure_dir_sizes())
    return out


_DIR_CACHE: Dict[str, Any] = {}


def _structure_dir_sizes() -> Dict[str, Any]:
    """File count and total size for the two structure directories.

    Stat'ing 1,100 files was the largest remaining cost of a poll once the queries were
    cached, and the answer only changes when a structure is written. The directory's own
    mtime moves on any create or delete, so it is enough of a key -- a file edited in place
    without changing the listing would be missed, which cannot happen here because these
    files are only ever written whole.
    """
    out: Dict[str, Any] = {}
    for label, path in (("structures", config.STATIC_DIR / "structures"),
                        ("reference_structures", config.STATIC_DIR / "reference_structures")):
        try:
            mtime = path.stat().st_mtime_ns
        except OSError:
            continue
        hit = _DIR_CACHE.get(label)
        if not hit or hit[0] != mtime:
            try:
                files = list(path.iterdir())
                hit = (mtime, len(files), round(sum(f.stat().st_size for f in files) / 1e6, 1))
            except OSError:
                continue
            _DIR_CACHE[label] = hit
        out[f"{label}_files"], out[f"{label}_mb"] = hit[1], hit[2]
    return out


def _db_signature() -> tuple:
    """What changes when, and only when, the pipeline writes.

    SQLite in WAL mode leaves the main file's mtime alone between checkpoints, so the WAL
    is stat'd too. Missing files count as zeros rather than raising: the signature only has
    to differ across a write, not describe the database.
    """
    sig = []
    for name in ("pants.db", "pants.db-wal"):
        try:
            st = (config.ROOT_DIR / name).stat()
            sig.append((st.st_mtime_ns, st.st_size))
        except OSError:
            sig.append((0, 0))
    return tuple(sig)


def _db_stats() -> Dict[str, Any]:
    marks = ",".join("?" * len(config.MEASURED_TIERS))
    tiers = list(config.MEASURED_TIERS)

    with connect() as c:
        # ---- sequence space searched -------------------------------------------------
        envs = _rows(c, """
            SELECT json_extract(params_json,'$.environment') env,
                   COALESCE(SUM(n_input),0) scanned, COUNT(*) runs
            FROM runs WHERE stage='recall' AND params_json IS NOT NULL
            GROUP BY 1 ORDER BY scanned DESC""")
        cand_env = {r["env"]: r["n"] for r in _rows(c, """
            SELECT source_environment env, COUNT(*) n FROM candidates GROUP BY 1""")}
        for e in envs:
            e["candidates"] = cand_env.get(e["env"], 0)
            e["per_million"] = (round(1e6 * e["candidates"] / e["scanned"], 1)
                                if e["scanned"] else None)

        # ---- structures --------------------------------------------------------------
        ref_src = {r["source"]: r["n"] for r in _rows(c, """
            SELECT source, COUNT(*) n FROM reference_structures GROUP BY 1""")}
        stats = {
            "sequences_scanned": _one(c, "SELECT COALESCE(SUM(n_input),0) FROM runs WHERE stage='recall'"),
            "environments": envs,

            "candidates": _one(c, "SELECT COUNT(*) FROM candidates"),
            "candidates_folded": _one(c, "SELECT COUNT(*) FROM structures"),
            "candidates_deferred": _one(c, "SELECT COUNT(*) FROM candidates WHERE structure_deferred=1"),
            "candidate_geometry": _one(c, "SELECT COUNT(*) FROM geometry WHERE triad_ser_resnum IS NOT NULL"),
            "candidate_plddt_mean": _one(c, "SELECT ROUND(AVG(plddt_mean),1) FROM structures"),
            "candidate_cleft_mean": _one(c, "SELECT ROUND(AVG(cleft_width_A),1) FROM geometry"),

            "reference_structures": sum(ref_src.values()),
            "reference_by_source": ref_src,
            "reference_with_triad": _one(c, "SELECT COUNT(*) FROM reference_geometry WHERE triad_ser_resnum IS NOT NULL"),
            "reference_knockouts": _one(c, "SELECT COUNT(*) FROM reference_structures WHERE catalytic_knockout_json IS NOT NULL"),
            "best_resolution": _one(c, "SELECT MIN(resolution_A) FROM reference_structures"),
            "mean_plddt_predicted": _one(c, "SELECT ROUND(AVG(plddt_mean),1) FROM reference_structures WHERE plddt_mean IS NOT NULL"),

            # ---- the reference set -----------------------------------------------------
            "characterised_total": _one(c, "SELECT COUNT(*) FROM characterised_enzymes"),
            "positives": _one(c, "SELECT COUNT(*) FROM characterised_enzymes WHERE is_positive=1"),
            "measured": _one(c, f"SELECT COUNT(*) FROM characterised_enzymes WHERE is_positive=1 AND source_ref IN ({marks})", tiers),
            "annotated_only": _one(c, f"SELECT COUNT(*) FROM characterised_enzymes WHERE is_positive=1 AND source_ref NOT IN ({marks})", tiers),
            "negatives": _one(c, "SELECT COUNT(*) FROM characterised_enzymes WHERE is_negative=1"),
            "near_misses": _one(c, "SELECT COUNT(*) FROM characterised_enzymes WHERE is_near_miss=1"),
            "within_family_negatives": _one(c, "SELECT COUNT(*) FROM characterised_enzymes WHERE source_ref='PAZy-nonPET'"),
            "named_enzymes": _one(c, "SELECT COUNT(*) FROM characterised_enzymes WHERE enzyme_id NOT LIKE '%:%' AND is_positive=1"),
            "tiers": _rows(c, """SELECT source_ref tier, COUNT(*) n FROM characterised_enzymes
                                 WHERE is_positive=1 GROUP BY 1 ORDER BY n DESC"""),

            # ---- evidence --------------------------------------------------------------
            "measurements": _one(c, "SELECT COUNT(*) FROM activity_measurements"),
            "measurements_cited": _one(c, "SELECT COUNT(*) FROM activity_measurements "
                                          "WHERE source_doi IS NOT NULL AND TRIM(source_doi) != ''"),
            "measurement_kinds": _rows(c, """SELECT parameter_type kind, COUNT(*) n
                                             FROM activity_measurements GROUP BY 1 ORDER BY n DESC"""),
            "evidence_codes": _rows(c, """SELECT evidence_code code, COUNT(*) n
                                          FROM activity_measurements WHERE evidence_code IS NOT NULL
                                          GROUP BY 1 ORDER BY n DESC"""),
            "distinct_papers": _one(c, "SELECT COUNT(DISTINCT source_doi) FROM activity_measurements "
                                       "WHERE source_doi IS NOT NULL AND TRIM(source_doi) != ''"),

            # ---- provenance ------------------------------------------------------------
            "stages": _rows(c, """SELECT stage, COUNT(*) runs, COALESCE(SUM(n_input),0) input,
                                         COALESCE(SUM(n_output),0) output, MAX(finished_at) last
                                  FROM runs GROUP BY 1 ORDER BY runs DESC"""),
            "runs_total": _one(c, "SELECT COUNT(*) FROM runs"),
            "wall_time_h": round((_one(c, "SELECT COALESCE(SUM(wall_time_s),0) FROM manifests") or 0) / 3600, 1),
            "manifests": _one(c, "SELECT COUNT(*) FROM manifests"),
            "recent_runs": _rows(c, """SELECT stage, label, status, started_at, finished_at,
                                              n_input, n_output FROM runs
                                       ORDER BY id DESC LIMIT 12"""),
            "sources": _rows(c, """SELECT name, version, retrieved_at, n_records, license, source_url
                                   FROM data_sources ORDER BY name"""),
            "training": _rows(c, """SELECT head_name, model_version, n_positives, n_negatives,
                                           auc, average_precision, brier_score,
                                           composition_baseline_auc, n_positive_clusters,
                                           evidence_level, trained_at, config_json
                                    FROM training_runs ORDER BY run_id DESC"""),
            "schema_version": SCHEMA_VERSION,
            "data_version": config.DATA_VERSION,
        }

        tv: Dict[str, str] = {}
        for r in _rows(c, "SELECT tool_versions_json FROM manifests WHERE tool_versions_json IS NOT NULL"):
            try:
                tv.update(json.loads(r["tool_versions_json"]) or {})
            except (TypeError, ValueError):
                pass
        stats["tools"] = dict(sorted(tv.items()))
        stats["models"] = [r["model_version"] for r in _rows(
            c, "SELECT DISTINCT model_version FROM manifests WHERE model_version IS NOT NULL")]
        stats["git_commit"] = _one(c, "SELECT git_commit FROM manifests WHERE git_commit IS NOT NULL ORDER BY id DESC LIMIT 1")

    return stats


_CACHE: Dict[str, Any] = {"signature": None, "stats": None}
_CACHE_LOCK = threading.Lock()


def gather() -> Dict[str, Any]:
    """The whole payload.

    Forty-odd queries cost about 17 ms, which is nothing once but is paid on every poll by
    every open tab, and between pipeline runs every one of those polls returns exactly what
    the last did. So the database half is recomputed only when the database has actually
    been written, keyed on its mtime and size; an idle poll then costs a couple of stat
    calls. Host metrics are read every time regardless -- they are a few /proc reads, and
    they are the half that genuinely does change between requests.

    Not a time-based cache on purpose: a TTL would either lag a run that just finished or
    expire for nothing while the pipeline is idle, and the file already says which it is.
    """
    sig = _db_signature()
    with _CACHE_LOCK:
        cached = _CACHE["stats"] if _CACHE["signature"] == sig else None
        if cached is None:
            cached = _db_stats()
            _CACHE["signature"], _CACHE["stats"] = sig, cached
    stats = dict(cached)
    stats["host"] = host_health()
    stats["generated_at"] = time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()) + " UTC"
    stats["generated_epoch"] = int(time.time())
    stats["data_changed_epoch"] = int(max(s[0] for s in sig) / 1e9)
    return stats


@bp.route("/api/stats")
def api_stats():
    """The whole payload, recomputed per request. The page polls this."""
    return jsonify(gather())


# Plain-language glosses. Kept beside the code that names these stages rather than in the
# template, so a stage renamed in the pipeline is renamed in one place.
STAGE_BLURB = {
    "recall": "search metagenome proteins for the polyesterase fold, then keep only those "
              "whose catalytic residues meet in space",
    "reference_structures": "fetch or fold a structure for each characterised enzyme and "
                            "superpose it onto IsPETase",
    "seeds": "load the curated wild types and derive each engineered variant from its "
             "parent plus a checked mutation list",
    "pazy": "import enzymes with published, measured activity on a plastic",
    "embed": "turn every sequence into 480 numbers with a frozen protein language model",
    "structure": "fold candidates with ESMFold and measure the active site",
    "positives_family": "add the wider enzyme family as similarity-labelled positives",
    "negatives": "collect hard negatives: the same fold, no PET activity",
    "harmonise": "reconcile activity values reported under different assays",
}

KIND_LABEL = {
    "topt": ("T<sub>opt</sub>", "the temperature at which the enzyme works fastest"),
    "km": ("K<sub>M</sub>", "how much substrate it takes to half-saturate the enzyme"),
    "kcat": ("k<sub>cat</sub>", "how many reactions one enzyme molecule runs per second"),
    "kcat_km": ("k<sub>cat</sub>/K<sub>M</sub>", "catalytic efficiency: the two above combined"),
    "tm": ("T<sub>m</sub>", "the temperature at which the fold falls apart"),
    "specific_activity": ("Specific activity", "product made per milligram of enzyme per minute"),
    "product_release": ("Product release", "how much PET breakdown product appeared"),
    "weight_loss": ("Weight loss", "how much of the plastic film disappeared"),
    "ph_opt": ("pH<sub>opt</sub>", "the acidity at which it works fastest"),
    "relative_activity": ("Relative activity", "activity as a percentage of a reference enzyme"),
    "performance_claim": ("Performance claim",
                          "a stated result from the paper — \"6-fold faster than the wild type\" — "
                          "kept as text because it has no single unit"),
    "catalytic_activity": ("Catalytic activity",
                           "product formed under the assay conditions the paper used"),
    "ordinal_activity": ("Ordinal activity",
                         "a rank rather than a number: this enzyme beat that one, with no value given"),
}


SOURCE_BLURB = {
    "UniProt": "curated protein sequences, annotations and catalytic-site positions",
    "PAZy": "the plastics-active enzyme database: an enzyme is listed because activity on "
            "a plastic was measured and published",
    "HGMP-SciDB": "human gut metagenome polyesterases from the deposit accompanying the paper",
}


def training_note(t: Dict[str, Any]) -> str:
    """One line saying what a training run actually showed, including the failures."""
    try:
        cfg = json.loads(t.get("config_json") or "{}")
    except (TypeError, ValueError):
        cfg = {}
    if t.get("auc") is None:
        return "not evaluable — too few independent clusters to split on"
    auc, base = t["auc"], t.get("composition_baseline_auc")
    if auc >= 0.999:
        return "perfect, and meaningless: it reproduced the similarity rule that made the labels"
    if "mixed" in (t.get("evidence_level") or ""):
        return "trained on similarity-derived labels"
    if auc < 0.6:
        return "at chance — the hardest contrast, and the honest answer"
    if base and auc - base < 0.1:
        return "barely clears amino-acid composition"
    return "clears the composition baseline by a real margin"


@bp.route("/stats")
def stats_page():
    return render_template("stats.html", active="stats", stats=gather(),
                           stage_blurb=STAGE_BLURB, source_blurb=SOURCE_BLURB,
                           training_note=training_note, kind_label=KIND_LABEL,
                           measured_tiers=set(config.MEASURED_TIERS),
                           env_label=ENV_LABEL)
