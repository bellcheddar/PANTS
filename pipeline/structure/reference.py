"""Structures and active-site geometry for the CHARACTERISED enzymes.

The candidate pipeline folds metagenomic hits and measures their active sites. This does
the same for the named reference enzymes, so the two can be put side by side: a candidate's
cleft width means very little on its own, and a great deal next to IsPETase's and LCC's
measured on the same code path.

**Coordinates come from three places, and which one matters.** An experimental entry is
used wherever one exists, because a crystal structure is evidence and a prediction is a
hypothesis. Failing that, AlphaFold, where the enzyme has a UniProt accession. The
engineered variants have neither -- DuraPETase, DepoPETase, LCC-A2, ThermoPETase and
FAST-PETase are mutation sets applied to a parent, with no entry of their own -- so those
are folded with the same ESMFold the candidates use. `source` records which, per enzyme,
because a geometric comparison across a mixture of crystal structures and predictions is
only honest if the reader can see which is which.

Everything is superposed onto IsPETase 6EQE at write time, exactly as the candidates are,
so the viewer can overlay any combination without aligning anything in the browser.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .. import config, http
from ..db import connect, now, retry_write
from ..db.manifest import stage_manifest
from . import fold, geometry

STAGE = "reference_structures"

# Served statically, alongside the candidate structures the viewer already reads.
REF_DIR = config.STATIC_DIR / "reference_structures"

PDB_COORD_URL = "https://files.rcsb.org/download/{pdb_id}.pdb"
ALPHAFOLD_API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"

# MHETase is in the reference set as the second enzyme of the two-step pathway, not as a
# polyesterase. It is Tannase-family, 603 aa, and its active site is not comparable to a
# PETase cleft, so it is carried for completeness and excluded from geometry comparisons.
NOT_A_POLYESTERASE = {"MHETase"}


def _fetch_pdb(pdb_id: str) -> Optional[str]:
    resp = http.get(PDB_COORD_URL.format(pdb_id=pdb_id.upper()))
    if resp is None or resp.status_code != 200 or not resp.text.startswith(("HEADER", "ATOM", "REMARK")):
        return None
    return resp.text


def _fetch_alphafold(acc: str) -> Optional[str]:
    """AlphaFold model by accession, resolving the file URL from the API.

    The URL is resolved rather than composed: building "AF-{acc}-F1-model_v4.pdb" by hand
    returned 404 for every entry once the database moved to v6.
    """
    meta = http.get_json(ALPHAFOLD_API.format(acc=acc))
    if not meta:
        return None
    url = (meta or [{}])[0].get("pdbUrl")
    if not url:
        return None
    resp = http.get(url)
    return resp.text if resp is not None and resp.status_code == 200 else None


def _mean_plddt(pdb_text: str) -> Optional[float]:
    """Mean CA pLDDT, always on the 0 to 100 convention.

    The two sources do not agree on scale and nothing upstream reconciles them. AlphaFold
    writes 0 to 100 into the B-factor column; ESMFold writes 0 to 1. `fold.fold_one` scales
    its RETURN value by PLDDT_SCALE but does not touch the PDB text, and this reads the
    text, so the raw numbers arriving here are 0-1 from one source and 0-100 from the
    other. Stored unreconciled, the column held 0.96 for an ESMFold model beside 92.01 for
    an AlphaFold one, and the page rendered "1.0" as a confidence score.

    Detected rather than assumed from the source argument: a B-factor column that never
    exceeds 1.0 is a fraction, since a real pLDDT of 1/100 would mean every residue was
    predicted as badly as it is possible to predict one.
    """
    vals = [float(l[60:66]) for l in pdb_text.splitlines()
            if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    if not vals:
        return None
    m = sum(vals) / len(vals)
    if max(vals) <= 1.0:
        m *= fold.PLDDT_SCALE
    return round(m, 2)


def _resolution(pdb_text: str) -> Optional[float]:
    """Resolution in Angstroms from the REMARK 2 record.

    Parsed from the text AFTER "RESOLUTION.", not from the first number on the line: the
    line reads `REMARK   2 RESOLUTION.    1.44 ANGSTROMS.` and the first float on it is the
    remark number 2. Taking that gave every experimental structure a resolution of exactly
    2.0 A, which is plausible enough to pass unnoticed.
    """
    for line in pdb_text.splitlines():
        if line.startswith("REMARK   2 RESOLUTION."):
            tail = line.split("RESOLUTION.", 1)[1]
            for tok in tail.split():
                try:
                    return float(tok)
                except ValueError:
                    continue
    return None


def sources_for(enzyme_id: str, uniprot: Optional[str], pdb_ids: List[str],
                sequence: Optional[str]) -> List[Tuple[str, str]]:
    """Ordered list of (source, id) to try. Experimental first, prediction last."""
    out: List[Tuple[str, str]] = []
    for pid in pdb_ids:
        out.append(("pdb", pid))
    if uniprot:
        out.append(("alphafold", uniprot))
    if sequence:
        out.append(("esmfold", enzyme_id))
    return out


def build(only: Optional[List[str]] = None, label: str = "v1") -> Dict[str, object]:
    """Fetch or fold, superpose onto 6EQE, measure, and write a viewer PDB per enzyme."""
    REF_DIR.mkdir(parents=True, exist_ok=True)
    report: Dict[str, object] = {"built": [], "skipped": [], "failed": []}

    with stage_manifest(STAGE, label=label) as m:
        with connect() as c:
            rows = c.execute(
                "SELECT enzyme_id, uniprot, pdb_ids_json, sequence, seq_length "
                "FROM characterised_enzymes "
                "WHERE enzyme_id NOT LIKE 'PAZy:%' AND enzyme_id NOT LIKE 'ESTHER:%' "
                "  AND enzyme_id NOT LIKE 'EC:%' AND enzyme_id NOT LIKE 'PLC:%' "
                "  AND is_positive=1 AND sequence IS NOT NULL "
                "ORDER BY enzyme_id").fetchall()

        targets = [r for r in rows if not only or r[0] in only]
        ref_cif = config.STATIC_DIR / "reference" / f"{config.ISPETASE_REFERENCE_PDB}.cif"

        for enzyme_id, uniprot, pdb_json, sequence, seq_len in targets:
            pdb_ids = json.loads(pdb_json or "[]")
            got = None
            for source, sid in sources_for(enzyme_id, uniprot, pdb_ids, sequence):
                if source == "pdb":
                    text = _fetch_pdb(sid)
                elif source == "alphafold":
                    text = _fetch_alphafold(sid)
                else:
                    text = _esmfold(sequence)
                if text:
                    got = (source, sid, text)
                    break

            if got is None:
                report["failed"].append(f"{enzyme_id}: no coordinates from any source")
                continue
            source, sid, pdb_text = got

            cif, rmsd, _frac = fold.superpose_onto_reference(pdb_text, ref_cif)
            if cif is None:
                # Superposition failing is informative rather than fatal: MHETase is a
                # different fold and does not align to a PETase.
                report["skipped"].append(f"{enzyme_id}: would not superpose onto "
                                         f"{config.ISPETASE_REFERENCE_PDB} ({source} {sid})")
                cif, rmsd = None, None

            dest = REF_DIR / f"{enzyme_id.replace('/', '_').replace('*', 's')}.pdb"
            written = (fold.write_viewer_pdb(cif, dest) if cif
                       else _write_plain(pdb_text, dest))
            if written is None:
                report["failed"].append(f"{enzyme_id}: could not write viewer coordinates")
                continue

            site = geometry.measure(dest)
            n_res = len({l[22:27] for l in dest.read_text().splitlines()
                         if l.startswith("ATOM")})
            offset = sequence_offset(sequence, site)
            if offset is None and site.ser_resnum is not None:
                report["skipped"].append(
                    f"{enzyme_id}: no single offset maps {source} {sid}'s numbering onto "
                    f"the stored sequence; the sequence panel cannot mark its triad")
            _persist(enzyme_id, source, sid, dest.name, pdb_text, rmsd, site, n_res,
                     offset)
            report["built"].append((enzyme_id, source, sid, n_res,
                                    site.cleft_width_A, site.triad_is_connected))

        m.counts(n_input=len(targets), n_output=len(report["built"]),
                 n_discarded=len(report["failed"]) + len(report["skipped"]))
    return report


def _write_plain(pdb_text: str, dest: Path) -> Optional[Path]:
    """Fallback for a structure that would not superpose: keep it, unaligned, and say so."""
    ss = fold.secondary_structure(pdb_text)
    dest.write_text(fold._ss_records(ss) + pdb_text)
    return dest


_MODEL = None


def _esmfold(sequence: str) -> Optional[str]:
    """Fold one sequence, loading ESMFold once per process.

    Loading dominates the cost (the model is 8.4 GB), so the handle is cached at module
    level: the reference set needs only a handful of folds and reloading per enzyme would
    take longer than the folding.
    """
    global _MODEL
    if len(sequence) > fold.MAX_FOLD_LENGTH:
        return None
    try:
        if _MODEL is None:
            _MODEL = fold.load_model()
        tok, model = _MODEL
        pdb_text, _plddt = fold.fold_one(tok, model, sequence)
        return pdb_text
    except Exception:
        return None


def sequence_offset(sequence: Optional[str], site: geometry.ActiveSite) -> Optional[int]:
    """How the structure's residue numbering relates to the stored sequence.

    `structure_resnum = sequence_position + offset`.

    A deposited construct is numbered by whoever deposited it, and that need not match a
    UniProt precursor: 7CEF is out by +42 against Cut190, 7QVH by +26 against HotPETase and
    4CG1 by -40 against TfCut2. Left unreconciled the viewer and the sequence panel
    disagree about which residues are catalytic, and the sequence panel is the one that is
    wrong, silently, on a page whose whole point is showing where the chemistry happens.

    Found by asking which single shift puts the MEASURED triad on a Ser, a His and an
    Asp/Glu in the stored sequence. Three independent positions agreeing pins it; a
    structure needing an indel rather than a shift returns None rather than a guess.
    """
    if not sequence or site.ser_resnum is None:
        return None
    for off in sorted(range(-120, 121), key=abs):
        ps, ph, pa = (site.ser_resnum - off, site.his_resnum - off, site.asp_resnum - off)
        if not all(0 < p <= len(sequence) for p in (ps, ph, pa)):
            continue
        if sequence[ps - 1] == "S" and sequence[ph - 1] == "H" and sequence[pa - 1] in "DE":
            return off
    return None


def _persist(enzyme_id: str, source: str, source_id: str, coord_name: str,
             raw_text: str, rmsd: Optional[float], site: geometry.ActiveSite,
             n_res: int, seq_offset: Optional[int] = None) -> None:
    def _do() -> None:
        with connect() as c:
            c.execute(
                "INSERT INTO reference_structures (enzyme_id, source, source_id, "
                " coord_path, plddt_mean, resolution_A, rmsd_ca_to_ispetase_A, "
                " superposition_reference, n_residues, built_at, model_version, "
                " seq_offset) VALUES (?,?,?,?,?,?,?,?,?,?,?,?) "
                "ON CONFLICT(enzyme_id) DO UPDATE SET source=excluded.source, "
                " source_id=excluded.source_id, coord_path=excluded.coord_path, "
                " plddt_mean=excluded.plddt_mean, resolution_A=excluded.resolution_A, "
                " rmsd_ca_to_ispetase_A=excluded.rmsd_ca_to_ispetase_A, "
                " n_residues=excluded.n_residues, built_at=excluded.built_at, "
                " seq_offset=excluded.seq_offset",
                (enzyme_id, source, source_id, coord_name,
                 _mean_plddt(raw_text) if source != "pdb" else None,
                 _resolution(raw_text) if source == "pdb" else None,
                 rmsd, config.ISPETASE_REFERENCE_PDB, n_res, now(),
                 config.ESMFOLD_MODEL if source == "esmfold" else source, seq_offset))
            c.execute(
                "INSERT INTO reference_geometry (enzyme_id, triad_ser_resnum, "
                " triad_his_resnum, triad_asp_resnum, ser_og_his_ne2_dist_A, "
                " his_nd1_asp_od_dist_A, ser_his_asp_angle_deg, oxyanion_n1_dist_A, "
                " oxyanion_n1_resnum, oxyanion_n2_dist_A, oxyanion_n2_resnum, "
                " oxyanion_n2_angle_deg, cleft_width_A, cleft_depth_A, n_cleft_residues, "
                " aromatic_clamp_residues_json) VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) "
                "ON CONFLICT(enzyme_id) DO UPDATE SET "
                " triad_ser_resnum=excluded.triad_ser_resnum, "
                " triad_his_resnum=excluded.triad_his_resnum, "
                " triad_asp_resnum=excluded.triad_asp_resnum, "
                " ser_og_his_ne2_dist_A=excluded.ser_og_his_ne2_dist_A, "
                " his_nd1_asp_od_dist_A=excluded.his_nd1_asp_od_dist_A, "
                " ser_his_asp_angle_deg=excluded.ser_his_asp_angle_deg, "
                " oxyanion_n1_dist_A=excluded.oxyanion_n1_dist_A, "
                " oxyanion_n1_resnum=excluded.oxyanion_n1_resnum, "
                " oxyanion_n2_dist_A=excluded.oxyanion_n2_dist_A, "
                " oxyanion_n2_resnum=excluded.oxyanion_n2_resnum, "
                " oxyanion_n2_angle_deg=excluded.oxyanion_n2_angle_deg, "
                " cleft_width_A=excluded.cleft_width_A, "
                " cleft_depth_A=excluded.cleft_depth_A, "
                " n_cleft_residues=excluded.n_cleft_residues, "
                " aromatic_clamp_residues_json=excluded.aromatic_clamp_residues_json",
                (enzyme_id, site.ser_resnum, site.his_resnum, site.asp_resnum,
                 site.ser_og_his_ne2_A, site.his_nd1_asp_od_A, site.ser_his_asp_angle_deg,
                 site.oxyanion_n1_A, site.oxyanion_n1_resnum, site.oxyanion_n2_A,
                 site.oxyanion_n2_resnum, site.oxyanion_n2_angle_deg,
                 site.cleft_width_A, site.cleft_depth_A, site.n_cleft_residues,
                 json.dumps(site.aromatic_clamp)))
    retry_write(_do)
