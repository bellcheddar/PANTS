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

# Shown by default in the Superpose view: the wild type everything is aligned onto, and
# one engineered variant that measurably works, so the frame opens with a comparison rather
# than a single structure.
DEFAULT_COMPARE = ["IsPETase", "FAST-PETase"]

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


def asset(filename: str) -> str:
    """Cache-busting version for a static file: its modification time.

    DATA_VERSION was doing this job and could not: it defaults to "dev" and never changes,
    so every asset URL was permanently `?v=dev` while nginx served /static/ as
    `max-age=31536000, immutable`. `immutable` is a promise to the browser that the body
    at this URL will never differ, so it does not revalidate at all -- for a year. Deploys
    landed correctly on the server and were invisible in the browser, with no error and
    nothing stale-looking to notice.

    Keyed on mtime, the URL changes exactly when the file does, which is what makes the
    immutable header safe rather than a trap.
    """
    path = config.STATIC_DIR / filename
    try:
        return str(int(path.stat().st_mtime))
    except OSError:
        return config.DATA_VERSION


@bp.app_context_processor
def inject_globals() -> Dict[str, Any]:
    with connect() as conn:
        try:
            n_struct = int(conn.execute("SELECT COUNT(*) FROM structures").fetchone()[0])
        except Exception:
            n_struct = 0
    return {"data_version": config.DATA_VERSION, "schema_version": SCHEMA_VERSION,
            "n_structures": n_struct, "asset": asset,
            "citation_links": citation_links}


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
    # The reference enzymes, ordered from the therapeutic end of the temperature range to
    # the industrial one, because that ordering IS the finding: the engineered lineage runs
    # away from 37 degC, and a table sorted by name would hide it.
    known = _named_enzyme_rows()
    for row in known:
        m = _mutations_for(row["enzyme_id"])
        row["n_mutations"] = (len(m["mutations"]) if m and m.get("mutations")
                              else (m or {}).get("n_substitutions"))
        row["parent"] = (m or {}).get("parent") or row.get("matched_positive_id")
        row["pdb_ids"] = json.loads(row.get("pdb_ids_json") or "[]")
        # A published performance figure where there is one, the curated description
        # otherwise, so no row is blank.
        row["headline_text"] = row.get("performance") or row.get("headline")
        row["reference"] = (m or {}).get("reference")
        row["reference_doi"] = _reference_doi(row["reference"])

    # Grouped by lineage, wild type first within each group, IsPETase's lineage first
    # overall. Sorting by name would scatter a lineage across the table; this keeps a wild
    # type next to the variants built on it, which is how the set is actually read.
    lineage_first = {"IsPETase": 0, "LCC": 1, "BhrPETase": 2, "Cut190": 3, "TfCut2": 4}
    # One colour per family, carried as a left rule on every row of that lineage. It
    # replaces a green bar that marked "Topt at or below 40 degC" and explained itself
    # nowhere: a coloured edge reads as a grouping, so it should encode the grouping.
    family_colour = {
        "IsPETase":  "#00d084",   # green
        "LCC":       "#4a9fd4",   # blue
        "BhrPETase": "#ff4d5e",   # red
        "Cut190":    "#fcb900",   # amber
        "TfCut2":    "#9b51e0",   # purple
        "MHETase":   "#6e78a6",   # grey: not a polyesterase
    }
    def _key(r):
        root = r.get("lineage_wt_id") or r["enzyme_id"]
        return (lineage_first.get(root, 9), root,
                0 if root == r["enzyme_id"] else 1,      # wild type heads its group
                -(r.get("identity_to_lineage_wt") or 0),  # then most-similar first
                r["enzyme_id"])
    known.sort(key=_key)

    # Grouped for the template: one block per lineage, the wild type first. The Parent
    # column is gone with this -- it said "wild type" for the wild types themselves, which
    # reads as though IsPETase has a parent called "wild type". Position in its own group
    # already says what a row is.
    groups: List[Dict[str, Any]] = []
    for row in known:
        root = row.get("lineage_wt_id") or row["enzyme_id"]
        if not groups or groups[-1]["root"] != root:
            groups.append({"root": root, "rows": []})
        row["is_root"] = (root == row["enzyme_id"])
        # HGMPs are five separate gut-metagenome enzymes, each its own lineage; one shared
        # colour keeps them legible as a set without inventing five more hues.
        row["family_colour"] = family_colour.get(
            root, "#ff6900" if root.startswith("HGMP") else "#3a4570")
        row["family_label"] = root
        groups[-1]["rows"].append(row)
    return render_template("home.html", active="home", counts=counts, known=known,
                           groups=groups,
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
    """The superposed viewer: reference enzymes and metagenomic candidates in one frame.

    Everything here was aligned onto IsPETase 6EQE when it was written, so any combination
    overlays directly and the browser aligns nothing. Two references are shown by default,
    IsPETase and FAST-PETase, because a candidate's geometry is only meaningful next to
    something known, and one reference cannot show how much a WORKING engineered enzyme
    already differs from the wild type.
    """
    with connect() as conn:
        cands = [dict(r) for r in conn.execute(
            "SELECT c.candidate_id AS sid, c.source_environment, c.seq_length, c.sequence, "
            "       s.plddt_mean, g.cleft_width_A, "
            "       g.triad_ser_resnum, g.triad_asp_resnum, g.triad_his_resnum "
            "FROM candidates c JOIN structures s ON s.candidate_id=c.candidate_id "
            "LEFT JOIN geometry g ON g.candidate_id=c.candidate_id "
            "ORDER BY s.plddt_mean DESC")]
        for r in cands:
            r["kind"] = "candidate"
            r["url"] = f"/static/structures/{r['sid']}.pdb"
            r["label"] = r["sid"].replace("PANTS-", "")

        refs = [dict(r) for r in conn.execute(
            "SELECT rs.enzyme_id AS sid, rs.coord_path, rs.source, rs.plddt_mean, "
            "       ce.sequence, ce.seq_length, ce.lineage_wt_id, "
            "       rg.cleft_width_A, rg.triad_ser_resnum, rg.triad_asp_resnum, "
            "       rg.triad_his_resnum "
            "FROM reference_structures rs "
            "JOIN characterised_enzymes ce ON ce.enzyme_id = rs.enzyme_id "
            "LEFT JOIN reference_geometry rg ON rg.enzyme_id = rs.enzyme_id "
            "WHERE rs.coord_path IS NOT NULL ORDER BY rs.enzyme_id")]
        # The two shown by default head the list: a picker whose checked entries are
        # scrolled out of sight looks like nothing is selected.
        refs.sort(key=lambda r: (DEFAULT_COMPARE.index(r["sid"])
                                 if r["sid"] in DEFAULT_COMPARE else len(DEFAULT_COMPARE),
                                 r["sid"]))
        for r in refs:
            r["kind"] = "reference"
            r["url"] = f"/static/reference_structures/{r['coord_path']}"
            r["label"] = r["sid"]
            m = _mutations_for(r["sid"])
            r["parent"] = (m or {}).get("parent")
            r["mutations"] = sorted(
                int("".join(ch for ch in x if ch.isdigit()))
                for x in (m or {}).get("mutations", [])
                if any(ch.isdigit() for ch in x))
            r["mut_labels"] = {int("".join(ch for ch in x if ch.isdigit())): x
                               for x in (m or {}).get("mutations", [])
                               if any(ch.isdigit() for ch in x)}

    for r in cands + refs:
        r["triad_at"] = {int(r[k]): lab for k, lab in
                         (("triad_ser_resnum", "Ser"), ("triad_his_resnum", "His"),
                          ("triad_asp_resnum", "Asp")) if r.get(k)}
        r.setdefault("mutations", [])
        r.setdefault("mut_labels", {})

    # ?id= accepts either kind, so a link from an enzyme page and a link from the
    # catalogue both land here with their subject already selected.
    requested = [x for x in request.args.getlist("id") if x]
    preselect = requested or DEFAULT_COMPARE
    return render_template("compare.html", active="compare", refs=refs, cands=cands,
                           preselect=preselect, default_refs=DEFAULT_COMPARE,
                           env_label=ENV_LABEL)


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


# ======================================================================================
# Known and characterised enzymes
#
# The reference set gets its own pages because it is the thing everything else is judged
# against: a candidate's cleft width or triad geometry means little on its own and a great
# deal beside IsPETase's, measured on the same code path. These pages are the one place
# where the sequence, the mutations, the structure, the measured activity and every
# external identifier for a named enzyme sit together.
# ======================================================================================

# Enzymes curated by name, as opposed to the bulk sets keyed by accession (PAZy:n,
# ESTHER:x, EC:x, PLC:x). Those are real data and are counted everywhere, but they are not
# what anyone means by "the known PETases".
# Qualified with the ce. alias: three joined tables carry an enzyme_id column and an
# unqualified reference is ambiguous.
NAMED_ENZYME_FILTER = (
    "ce.enzyme_id NOT LIKE 'PAZy:%' AND ce.enzyme_id NOT LIKE 'ESTHER:%' "
    "AND ce.enzyme_id NOT LIKE 'EC:%' AND ce.enzyme_id NOT LIKE 'PLC:%'"
)

# Experimental evidence outranks a review. Both are kept (UniProt curates IsPETase at
# 40 degC, the review gives 30-35, measured on different substrates), so the display has to
# choose deterministically rather than by whichever row the database returns first.
_EVIDENCE_RANK = "CASE a.evidence_code WHEN 'ECO:0000269' THEN 0 ELSE 1 END"


def _named_enzyme_rows() -> List[Dict[str, Any]]:
    with connect() as conn:
        # Wrapped in a subselect so the ORDER BY sorts the RESULT rather than
        # re-evaluating the correlated subqueries in the sort context, which produced two
        # separately-sorted groups: four enzymes with a Topt, then the NULLs, then the
        # remaining nine with a Topt sorted again from the start.
        return [dict(r) for r in conn.execute(f"""SELECT * FROM (
            SELECT ce.enzyme_id, ce.uniprot, ce.organism, ce.seq_length, ce.source_ref,
                   ce.pdb_ids_json, ce.matched_positive_id, ce.family,
                   ce.activity_substrate_notes, ce.headline,
                   ce.lineage_wt_id, ce.identity_to_lineage_wt,
                   (SELECT a.rate_value FROM activity_measurements a
                     WHERE a.enzyme_id=ce.enzyme_id AND a.parameter_type='topt'
                     ORDER BY {_EVIDENCE_RANK} LIMIT 1)               AS topt_c,
                   (SELECT a.raw_text FROM activity_measurements a
                     WHERE a.enzyme_id=ce.enzyme_id AND a.parameter_type='topt'
                     ORDER BY {_EVIDENCE_RANK} LIMIT 1)               AS topt_text,
                   (SELECT a.rate_value FROM activity_measurements a
                     WHERE a.enzyme_id=ce.enzyme_id AND a.parameter_type='ph_opt'
                     ORDER BY {_EVIDENCE_RANK} LIMIT 1)               AS ph_opt,
                   (SELECT a.raw_text FROM activity_measurements a
                     WHERE a.enzyme_id=ce.enzyme_id
                       AND a.parameter_type='performance_claim' LIMIT 1) AS performance,
                   rs.source AS struct_source, rs.source_id AS struct_source_id,
                   rs.coord_path, rs.plddt_mean, rs.resolution_A,
                   rs.rmsd_ca_to_ispetase_A, rs.n_residues,
                   rg.cleft_width_A, rg.cleft_depth_A, rg.n_cleft_residues,
                   rg.triad_ser_resnum, rg.triad_asp_resnum, rg.triad_his_resnum,
                   rg.ser_og_his_ne2_dist_A, rg.his_nd1_asp_od_dist_A,
                   rg.oxyanion_n1_dist_A, rg.oxyanion_n1_resnum,
                   rg.oxyanion_n2_dist_A, rg.oxyanion_n2_resnum,
                   rg.aromatic_clamp_residues_json
            FROM characterised_enzymes ce
            LEFT JOIN reference_structures rs ON rs.enzyme_id = ce.enzyme_id
            LEFT JOIN reference_geometry  rg ON rg.enzyme_id = ce.enzyme_id
            WHERE {NAMED_ENZYME_FILTER} AND ce.is_positive = 1
            ) ORDER BY (topt_c IS NULL), topt_c, enzyme_id""")]


def citation_links(source: Optional[str]) -> List[Dict[str, str]]:
    """Split a source field into linkable citations.

    The column holds several shapes, sometimes mixed in one value:
    `PMID:22194294;PMID:32269349`, `PMID:39551294;doi:10.1016/j.ijbiomac.2024.137732`,
    and bare URLs. Each token is resolved to its own link so a row citing three papers
    offers three links rather than one opaque string.

    Anything unrecognised is returned as text with no url, because a citation that cannot
    be resolved should look unresolved rather than link somewhere plausible and wrong.
    """
    out: List[Dict[str, str]] = []
    if not source:
        return out
    for tok in (x.strip() for x in source.replace(",", ";").split(";")):
        if not tok:
            continue
        low = tok.lower()
        if low.startswith("pmid:"):
            pmid = tok.split(":", 1)[1].strip()
            out.append({"label": f"PMID {pmid}",
                        "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"})
        elif low.startswith("doi:"):
            doi = tok.split(":", 1)[1].strip()
            out.append({"label": doi, "url": f"https://doi.org/{doi}"})
        elif low.startswith("http"):
            out.append({"label": "source", "url": tok})
        else:
            out.append({"label": tok, "url": ""})
    return out


def _reference_doi(reference: Optional[str]) -> Optional[str]:
    """Verified DOI for a citation string, or None.

    None is a real answer here, not a gap to paper over: the 2025 review and the Cut190
    reference have no DOI recorded, and a plausible-looking guess would send a reader to
    the wrong paper. Only DOIs checked through Crossref -- resolved, and the type, title,
    journal and year matched against the citation -- are in the map.
    """
    if not reference:
        return None
    try:
        from pipeline.recall import seeds
    except Exception:
        return None
    return seeds.REFERENCE_DOI.get(reference)


def _mutations_for(enzyme_id: str) -> Optional[Dict[str, Any]]:
    """Mutation set and parent, read from the curated seed definitions.

    Imported lazily and defensively: pipeline.recall.seeds is pure stdlib today, but the
    web venv is deliberately thin and this page must not be the thing that drags a
    dependency onto the droplet.
    """
    try:
        from pipeline.recall import seeds
    except Exception:
        return None
    # PDB-derived first: HotPETase and Cut190**SS appear in BOTH lists, as a deliberately
    # unconfirmed VARIANTS entry (no mutation list) and as a PDB_DERIVED entry carrying the
    # published substitution count. Checking VARIANTS first returned the empty list and the
    # count never surfaced, so the table showed a dash for an enzyme with 21 known
    # substitutions. The two are merged instead of one shadowing the other.
    pdb_derived = getattr(seeds, "PDB_DERIVED", {}).get(enzyme_id)

    for v in seeds.VARIANTS:
        if v.enzyme_id == enzyme_id:
            out = {"parent": v.parent, "mutations": v.mutations,
                   "confirmed": v.mutations_confirmed, "reference": v.reference,
                   "notes": v.notes}
            if pdb_derived and not v.mutations:
                pdb_id, parent, expected, ref = pdb_derived
                out.update({"n_substitutions": expected, "confirmed": True,
                            "reference": ref, "parent": parent,
                            "notes": f"Sequence taken from the deposited construct "
                                     f"{pdb_id}: {expected} substitutions against "
                                     f"{parent}."})
            return out
    for name, (pdb_id, parent, expected, ref) in getattr(seeds, "PDB_DERIVED", {}).items():
        if name == enzyme_id:
            # No mutation LIST -- the sequence came from a deposited construct -- but the
            # published substitution count is known and checked, so the table can show it.
            return {"parent": parent, "mutations": [], "n_substitutions": expected,
                    "confirmed": True, "reference": ref,
                    "notes": f"Sequence taken from the deposited construct {pdb_id}: "
                             f"{expected} substitutions against {parent}."}
    return None


def _links_for(row: Dict[str, Any], mut: Optional[Dict[str, Any]]) -> Dict[str, List[Dict[str, str]]]:
    """Every external identifier and cross-reference for one enzyme, grouped.

    Gathered here rather than in the template so a missing accession produces no link
    instead of a dead one pointing at a URL built from None.
    """
    pdb_ids = json.loads(row.get("pdb_ids_json") or "[]")
    out: Dict[str, List[Dict[str, str]]] = {"structure": [], "sequence": [],
                                            "activity": [], "internal": []}
    for pid in pdb_ids:
        out["structure"].append({"label": f"PDB {pid}",
                                 "url": f"https://www.rcsb.org/structure/{pid}"})
    src, sid = row.get("struct_source"), row.get("struct_source_id")
    if src == "alphafold" and sid:
        out["structure"].append({"label": f"AlphaFold {sid}",
                                 "url": f"https://alphafold.ebi.ac.uk/entry/{sid}"})
    if row.get("uniprot"):
        acc = row["uniprot"]
        out["sequence"].append({"label": f"UniProt {acc}",
                                "url": f"https://www.uniprot.org/uniprotkb/{acc}"})
        out["sequence"].append({"label": "InterPro domains",
                                "url": f"https://www.ebi.ac.uk/interpro/protein/UniProtKB/{acc}/"})
    # PAZy indexes by its own id, and the same protein is often present under an accession
    # we already hold, so the cross-reference is by accession rather than by name.
    if row.get("uniprot"):
        out["activity"].append({"label": "PAZy (search by accession)",
                                "url": f"https://pazy.eu/?s={row['uniprot']}"})
    out["activity"].append({"label": "ESTHER family database",
                            "url": "https://bioweb.supagro.inrae.fr/ESTHER/"})
    if mut and mut.get("parent"):
        out["internal"].append({"label": f"Parent enzyme: {mut['parent']}",
                                "url": f"/enzyme/{mut['parent']}"})
    return out


@bp.route("/enzyme/<path:enzyme_id>")
def enzyme(enzyme_id: str):
    """One characterised enzyme: everything known about it, in one place."""
    rows = {r["enzyme_id"]: r for r in _named_enzyme_rows()}
    row = rows.get(enzyme_id)
    if row is None:
        abort(404)

    with connect() as conn:
        seq = conn.execute("SELECT sequence FROM characterised_enzymes WHERE enzyme_id=?",
                           (enzyme_id,)).fetchone()
        measurements = [dict(r) for r in conn.execute(
            "SELECT parameter_type, rate_value, rate_units, temperature_c, ph, "
            "       substrate_form, raw_text, evidence_code, source_doi, "
            "       extraction_confidence "
            "FROM activity_measurements WHERE enzyme_id=? "
            "ORDER BY parameter_type, rate_value", (enzyme_id,))]
        # Derived enzymes that name this one as their parent: the lineage, downward.
        children = [r[0] for r in conn.execute(
            "SELECT enzyme_id FROM characterised_enzymes WHERE matched_positive_id=? "
            "ORDER BY enzyme_id", (enzyme_id,))]
        # Candidates whose nearest characterised enzyme is this one. This is the join that
        # makes the page a hub rather than a leaf: it connects the reference set to the
        # metagenomic catalogue.
        nearest = [dict(r) for r in conn.execute(
            "SELECT candidate_id, source_environment, seq_length, recall_bitscore, "
            "       recall_profile_identity "
            "FROM candidates WHERE nearest_characterised_id=? "
            "ORDER BY recall_bitscore DESC LIMIT 25", (enzyme_id,))]
        n_nearest = int(conn.execute(
            "SELECT COUNT(*) FROM candidates WHERE nearest_characterised_id=?",
            (enzyme_id,)).fetchone()[0])
        # Everything else with a structure, for the overlay picker.
        others = [dict(r) for r in conn.execute(
            "SELECT rs.enzyme_id, rs.coord_path, rs.source, rg.triad_ser_resnum, "
            "       rg.triad_asp_resnum, rg.triad_his_resnum "
            "FROM reference_structures rs "
            "LEFT JOIN reference_geometry rg ON rg.enzyme_id=rs.enzyme_id "
            "WHERE rs.enzyme_id != ? AND rs.coord_path IS NOT NULL "
            "ORDER BY rs.enzyme_id", (enzyme_id,))]

    mut = _mutations_for(enzyme_id)
    row = dict(row)
    row["reference_doi"] = _reference_doi((mut or {}).get("reference"))

    # Map 1-based sequence position -> mutation label, for marking up the sequence panel.
    # Every confirmed variant in this set matched its parent at offset 0, which is what
    # apply_mutations verifies, so the stated position indexes the stored sequence
    # directly. If a set ever needed an offset that fact would be recorded with it, and
    # this would have to apply it rather than assume zero.
    mut_at: Dict[int, str] = {}
    if mut and mut.get("mutations"):
        for label in mut["mutations"]:
            digits = "".join(ch for ch in label if ch.isdigit())
            if digits:
                mut_at[int(digits)] = label

    # Triad positions, so the sequence can mark them in the viewer's own yellow.
    triad_at: Dict[int, str] = {}
    for key, label in (("triad_ser_resnum", "Ser"), ("triad_his_resnum", "His"),
                       ("triad_asp_resnum", "Asp")):
        if row.get(key):
            triad_at[int(row[key])] = label

    return render_template(
        "enzyme.html", active="home", e=row, seq=seq[0] if seq else None,
        triad_at=triad_at,
        mut=mut, measurements=measurements, children=children, nearest=nearest,
        n_nearest=n_nearest, others=others, links=_links_for(row, mut),
        mut_at=mut_at,
        clamp=json.loads(row.get("aromatic_clamp_residues_json") or "[]"),
        env_label=ENV_LABEL)
