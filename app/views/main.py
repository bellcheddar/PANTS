"""The v1 tabs: Home, Catalogue, Candidate, Superpose, Methods.

Reads precomputed SQLite and serves static mmCIF. Nothing here imports torch: the droplet
has 3.8 GB shared with five other applications, and keeping this process thin is the whole
reason the pipeline and web virtual environments are separate.
"""

from __future__ import annotations

import json
from typing import Any, Dict, List, Optional

from flask import Blueprint, abort, jsonify, render_template, request

from pipeline import config
from pipeline.db import connect
from pipeline.db.schema import SCHEMA_VERSION

bp = Blueprint("main", __name__)

ENV_LABEL = {
    "compost": "Compost",
    "marine_plastisphere": "Marine plastisphere",
    "landfill": "Landfill",
    "wastewater": "Wastewater",
    "human_gut": "Human gut",
    "unknown": "Unknown",
}


# Tiers whose label rests on a measurement. PAZy is here because its inclusion criterion
# is that activity was measured on a plastic and published; the EC-auto-annotated bulk is
# NOT here because those labels were assigned by sequence similarity.
# Single definition, in pipeline.config: the app and the release bundle reported
# different measured counts once because each kept its own copy.
MEASURED_TIERS = config.MEASURED_TIERS
MEASURED_MARKS = ",".join("?" * len(MEASURED_TIERS))


def _counts() -> Dict[str, int]:
    with connect() as conn:
        def n(sql: str, params=()) -> int:
            try:
                return int(conn.execute(sql, params).fetchone()[0])
            except Exception:
                return 0
        return {
            "candidates": n("SELECT COUNT(*) FROM candidates"),
            "structures": n("SELECT COUNT(*) FROM structures"),
            "positives": n("SELECT COUNT(*) FROM characterised_enzymes WHERE is_positive=1"),
            # Every tier whose label rests on an EXPERIMENT rather than on similarity.
            # Defined once here and imported by anything that needs it, because deriving
            # this list twice is how the app came to report 12 where the README said 17.
            "evidenced": n("SELECT COUNT(*) FROM characterised_enzymes WHERE is_positive=1 "
                           f"AND source_ref IN ({MEASURED_MARKS})", MEASURED_TIERS),
            "negatives": n("SELECT COUNT(*) FROM characterised_enzymes WHERE is_negative=1"),
            "near_misses": n("SELECT COUNT(*) FROM characterised_enzymes WHERE is_near_miss=1"),
            "measurements": n("SELECT COUNT(*) FROM activity_measurements"),
            # Sum of what recall actually scanned. Stored per run rather than derived,
            # so it stays right when a collection is added or rescanned.
            "sequences": n("SELECT COALESCE(SUM(n_input),0) FROM runs WHERE stage='recall'"),
        }


@bp.app_context_processor
def inject_globals() -> Dict[str, Any]:
    with connect() as conn:
        try:
            n_struct = int(conn.execute("SELECT COUNT(*) FROM structures").fetchone()[0])
        except Exception:
            n_struct = 0
    return {"data_version": config.DATA_VERSION, "schema_version": SCHEMA_VERSION,
            "n_structures": n_struct}


@bp.route("/")
def home():
    counts = _counts()
    with connect() as conn:
        # The environment a recall run covered is in its params_json, so the scanned
        # totals join back on that rather than being recounted from the FASTA files.
        scanned = {r["env"]: r["n"] for r in conn.execute(
            "SELECT json_extract(params_json,'$.environment') env, COALESCE(SUM(n_input),0) n "
            "FROM runs WHERE stage='recall' AND params_json IS NOT NULL GROUP BY 1")}
        by_env = [dict(r) for r in conn.execute(
            "SELECT source_environment env, COUNT(*) n, "
            "       AVG(recall_profile_identity) mean_ident "
            "FROM candidates GROUP BY 1 ORDER BY n DESC")]
        for row in by_env:
            row["scanned"] = scanned.get(row["env"], 0)
            # Yield per million is the comparable figure: raw candidate counts compare
            # environments that were sampled to very different depths.
            row["per_million"] = (1e6 * row["n"] / row["scanned"]) if row["scanned"] else None
        top = [dict(r) for r in conn.execute(
            "SELECT c.candidate_id, c.source_environment, c.seq_length, "
            "       c.recall_bitscore, c.recall_profile_identity, c.nearest_characterised_id, "
            "       s.plddt_mean, g.cleft_width_A "
            "FROM candidates c "
            "LEFT JOIN structures s ON s.candidate_id=c.candidate_id "
            "LEFT JOIN geometry g ON g.candidate_id=c.candidate_id "
            "ORDER BY c.recall_bitscore DESC")]
    return render_template("home.html", active="home", counts=counts,
                           by_env=by_env, top=top, n_visible=10, env_label=ENV_LABEL)


@bp.route("/catalogue")
def catalogue():
    env = request.args.get("env") or ""
    with connect() as conn:
        sql = ("SELECT c.candidate_id, c.source_environment, c.seq_length, "
               "       c.recall_evalue, c.recall_bitscore, c.recall_profile_identity, "
               "       c.nearest_characterised_id, c.structure_deferred, s.plddt_mean, "
               "       g.cleft_width_A, g.triad_ser_resnum, g.triad_asp_resnum, g.triad_his_resnum "
               "FROM candidates c "
               "LEFT JOIN structures s ON s.candidate_id=c.candidate_id "
               "LEFT JOIN geometry g ON g.candidate_id=c.candidate_id ")
        params: List[Any] = []
        if env:
            sql += "WHERE c.source_environment=? "
            params.append(env)
        sql += "ORDER BY c.recall_bitscore DESC"
        rows = [dict(r) for r in conn.execute(sql, params)]
        envs = [r[0] for r in conn.execute(
            "SELECT DISTINCT source_environment FROM candidates ORDER BY 1")]
    return render_template("catalogue.html", active="catalogue", rows=rows,
                           envs=envs, env=env, env_label=ENV_LABEL)


@bp.route("/candidate/<cid>")
def candidate(cid: str):
    with connect() as conn:
        row = conn.execute(
            "SELECT c.*, s.plddt_mean, s.mmcif_path, s.rmsd_ca_to_ispetase_A, "
            "       s.tm_score_to_ispetase, s.structure_method "
            "FROM candidates c LEFT JOIN structures s ON s.candidate_id=c.candidate_id "
            "WHERE c.candidate_id=?", (cid,)).fetchone()
        if row is None:
            abort(404)
        geom = conn.execute("SELECT * FROM geometry WHERE candidate_id=?", (cid,)).fetchone()
    return render_template("candidate.html", active="catalogue", c=dict(row),
                           g=dict(geom) if geom else None, env_label=ENV_LABEL)


@bp.route("/compare")
def compare():
    """The superposed viewer. Structures are already in a common frame."""
    with connect() as conn:
        rows = [dict(r) for r in conn.execute(
            "SELECT c.candidate_id, c.source_environment, c.seq_length, "
            "       s.plddt_mean, g.cleft_width_A, "
            "       g.triad_ser_resnum, g.triad_asp_resnum, g.triad_his_resnum "
            "FROM candidates c JOIN structures s ON s.candidate_id=c.candidate_id "
            "LEFT JOIN geometry g ON g.candidate_id=c.candidate_id "
            "ORDER BY s.plddt_mean DESC")]
    preselect = [r for r in request.args.getlist("id") if r]
    return render_template("compare.html", active="compare", rows=rows,
                           preselect=preselect, env_label=ENV_LABEL)


@bp.route("/api/structures")
def api_structures():
    """Structure metadata for the viewer, including triad residue numbers."""
    ids = request.args.getlist("id")
    if not ids:
        return jsonify([])
    marks = ",".join("?" * len(ids))
    with connect() as conn:
        rows = [dict(r) for r in conn.execute(
            f"SELECT c.candidate_id, c.source_environment, s.mmcif_path, s.plddt_mean, "
            f"       g.triad_ser_resnum ser, g.triad_asp_resnum asp, g.triad_his_resnum his "
            f"FROM candidates c JOIN structures s ON s.candidate_id=c.candidate_id "
            f"LEFT JOIN geometry g ON g.candidate_id=c.candidate_id "
            f"WHERE c.candidate_id IN ({marks})", ids)]
    return jsonify(rows)


@bp.route("/methods")
def methods():
    """Full transparency page. Spec section 7 calls this non-negotiable given how easily
    this kind of tool is over-read, and this project has an unusual amount to disclose."""
    with connect() as conn:
        runs = [dict(r) for r in conn.execute(
            "SELECT stage, label, status, started_at, n_input, n_output, n_discarded "
            "FROM runs ORDER BY id DESC LIMIT 20")]
        tiers = [dict(r) for r in conn.execute(
            "SELECT source_ref tier, COUNT(*) n FROM characterised_enzymes "
            "WHERE is_positive=1 GROUP BY 1 ORDER BY n DESC")]
        training = [dict(r) for r in conn.execute(
            "SELECT * FROM training_runs ORDER BY run_id DESC")]
        sources = [dict(r) for r in conn.execute("SELECT * FROM data_sources ORDER BY name")]
        manifests = [dict(r) for r in conn.execute(
            "SELECT stage, model_version, schema_version, git_commit, wall_time_s, written_at "
            "FROM manifests ORDER BY id DESC LIMIT 12")]
    return render_template("methods.html", active="methods", runs=runs, tiers=tiers,
                           training=training, sources=sources, manifests=manifests,
                           counts=_counts())


@bp.route("/healthz")
def healthz():
    return {"status": "ok", "data_version": config.DATA_VERSION}
