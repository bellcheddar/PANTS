"""ESMFold prediction, superposition and geometry, written per candidate as it completes.

Three things this stage does deliberately.

**Resumable.** The model takes ~31 minutes to load and each structure ~110 seconds, so a
full run is hours. Every candidate is written to disk the moment it finishes and skipped on
a later pass, so an interrupted run resumes instead of restarting. A restart that begins
from zero turns a supervisor into a machine for losing work.

**Superposed at write time.** Every structure is aligned onto IsPETase (6EQE) before being
saved, so the mmCIF files ship already in a common frame. The web app can then load several
at once and have them overlay correctly without doing any alignment in the browser, which
is what makes an interactive superposed view cheap on the client.

**pLDDT rescaled to 0 to 100.** ESMFold's output here is on a 0 to 1 scale. The smoke test
printed "mean pLDDT 1.0", which reads as perfect confidence and is really 100 on the
conventional scale. Storing the raw 0 to 1 value would put a number in the database that
every convention says means something else.
"""

from __future__ import annotations

import json
import threading
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np

from .. import config
from ..db import connect, now, retry_write
from . import geometry

# TM-align's C bindings use static buffers and are not thread safe. AlphaFraud guards this
# with the same pattern; reuse it rather than rediscovering the corruption.
_TMALIGN_LOCK = threading.Lock()

PLDDT_SCALE = 100.0     # ESMFold emits 0 to 1 here; the convention is 0 to 100.

# Fold time is quadratic in length: measured seconds = 9.63e-04 * L^2.04 over 24 real
# folds on this machine. That makes the long tail dominate: 14 of 128 candidates exceed
# 450 aa and would consume 46% of the total compute.
#
# The scientific reason matters more than the scheduling one. A single-domain alpha/beta
# hydrolase is roughly 250 to 350 aa, so an 859-residue "candidate" is almost certainly a
# fusion, a multi-domain protein or a misassembly. The cleft measurement also assumes one
# active site in one domain, so its geometry would be misleading rather than merely slow.
#
# This is the SAME window already applied to the training set in pipeline/train/dataset.py.
# Applying a length sanity check to one side of the pipeline and not the other was an
# inconsistency, not a deliberate choice.
MAX_FOLD_LENGTH = 450


@dataclass
class FoldResult:
    candidate_id: str
    mmcif_path: str
    plddt_mean: float
    tm_score: Optional[float]
    rmsd_A: Optional[float]
    seconds: float


def load_model():
    """ESMFold on CPU. Loading dominates the run, so it is done once per process."""
    import torch
    from transformers import AutoTokenizer, EsmForProteinFolding

    tok = AutoTokenizer.from_pretrained(config.ESMFOLD_MODEL)
    model = EsmForProteinFolding.from_pretrained(config.ESMFOLD_MODEL, low_cpu_mem_usage=True)
    model.eval()
    torch.set_grad_enabled(False)
    if hasattr(model, "trunk") and hasattr(model.trunk, "set_chunk_size"):
        model.trunk.set_chunk_size(64)   # trades speed for peak memory in the trunk
    return tok, model


def fold_one(tok, model, sequence: str) -> Tuple[str, float]:
    """Returns (PDB text, mean pLDDT on a 0 to 100 scale)."""
    import torch

    with torch.no_grad():
        enc = tok([sequence], return_tensors="pt", add_special_tokens=False)
        out = model(**enc)
    pdb = model.output_to_pdb(out)[0]
    plddt = float(out["plddt"][0, :, 1].mean()) * PLDDT_SCALE
    return pdb, plddt


def superpose_onto_reference(pdb_text: str, reference_cif: Path
                             ) -> Tuple[Optional[str], Optional[float], Optional[float]]:
    """Align a predicted structure onto the reference and return superposed mmCIF.

    Returns (mmcif text, rmsd, tm_score_proxy). Superposing here rather than in the browser
    is what lets the Compare view overlay several candidates without client-side alignment.
    """
    import gemmi

    st = gemmi.read_structure_string(pdb_text, format=gemmi.CoorFormat.Pdb)
    st.setup_entities()
    ref = gemmi.read_structure(str(reference_cif))
    ref.setup_entities()
    ref.remove_ligands_and_waters()

    try:
        pol_s = st[0][0].get_polymer()
        pol_r = ref[0][0].get_polymer()
        sup = gemmi.calculate_superposition(pol_r, pol_s, pol_s.check_polymer_type(),
                                            gemmi.SupSelect.CaP)
        st[0].transform_pos_and_adp(sup.transform)
        rmsd = float(sup.rmsd)
        n_aligned = int(sup.count)
    except Exception:
        return None, None, None

    st.setup_entities()
    cif = st.make_mmcif_document().as_string()

    # A length-normalised alignment fraction. Not a true TM-score (that needs TM-align's
    # own optimisation), so it is stored under its own name rather than passed off as one.
    ref_len = max(1, sum(1 for _ in pol_r))
    frac = n_aligned / ref_len
    return cif, rmsd, frac


def secondary_structure(pdb_text: str) -> Dict[int, str]:
    """Per-residue secondary structure, as {residue number: 'h'|'s'|'c'}.

    ESMFold emits no HELIX or SHEET records, and neither gemmi's mmCIF nor a bare PDB then
    gives a viewer anything to build a cartoon from: 3Dmol renders a featureless tube and
    Mol* falls back to lines. biotite's P-SEA implementation assigns it from coordinates,
    which is the same route AlphaFraud takes for its ribbon viewer.
    """
    import io
    import biotite.structure as struc
    import biotite.structure.io.pdb as pdb_io

    try:
        f = pdb_io.PDBFile.read(io.StringIO(pdb_text))
        arr = f.get_structure(model=1)
        arr = arr[struc.filter_amino_acids(arr)]
        sse = struc.annotate_sse(arr)              # 'a' | 'b' | 'c' per residue
        res_ids = struc.get_residues(arr)[0]
    except Exception:
        return {}

    mapping: Dict[int, str] = {}
    for rid, code in zip(res_ids, sse):
        mapping[int(rid)] = {"a": "h", "b": "s"}.get(str(code), "c")
    return mapping


def _ss_records(ss: Dict[int, str]) -> str:
    """HELIX/SHEET records from a per-residue map, in real PDB columns.

    3Dmol's build does not set atom.ss from these itself, but the client parses them to do
    so (the same trick AlphaFraud uses). Writing them is what turns a tube into a cartoon.
    """
    if not ss:
        return ""
    lines: List[str] = []
    runs: List[Tuple[str, int, int]] = []
    prev_code, start, prev_num = None, None, None
    for num in sorted(ss):
        code = ss[num]
        if code != prev_code or (prev_num is not None and num != prev_num + 1):
            if prev_code in ("h", "s") and start is not None:
                runs.append((prev_code, start, prev_num))
            start, prev_code = num, code
        prev_num = num
    if prev_code in ("h", "s") and start is not None:
        runs.append((prev_code, start, prev_num))

    h = e = 0
    for code, a, b in runs:
        if b - a < 2:                       # a two-residue "helix" is noise
            continue
        if code == "h":
            h += 1
            lines.append(f"HELIX  {h:>3} {h:>3} ALA A {a:>4}  ALA A {b:>4}  1"
                         f"{'':30}{b-a+1:>5}")
        else:
            e += 1
            lines.append(f"SHEET  {e:>3} {e:>3} 1 ALA A {a:>4}  ALA A {b:>4}  0")
    return "\n".join(lines) + ("\n" if lines else "")


def already_done(candidate_id: str, out_dir: Path) -> bool:
    p = out_dir / f"{candidate_id}.cif"
    return p.exists() and p.stat().st_size > 0


def defer_long_candidates(max_length: int = MAX_FOLD_LENGTH) -> int:
    """Mark over-length candidates as deferred, with the reason on the row."""
    def _do() -> int:
        with connect() as conn:
            cur = conn.execute(
                "UPDATE candidates SET structure_deferred=1, structure_deferred_reason=? "
                "WHERE seq_length > ?",
                (f"length > {max_length} aa: probable fusion, multi-domain protein or "
                 f"misassembly. Fold time is quadratic in length and the cleft measurement "
                 f"assumes a single active site in one domain.", max_length))
            return cur.rowcount
    return retry_write(_do)


def run(candidates: Sequence[Tuple[str, str]], out_dir: Optional[Path] = None,
        reference_cif: Optional[Path] = None, progress_path: Optional[Path] = None,
        limit: Optional[int] = None, max_length: Optional[int] = MAX_FOLD_LENGTH
        ) -> List[FoldResult]:
    """Fold, superpose, measure and persist each candidate. Safe to re-run."""
    out_dir = out_dir or config.STRUCTURE_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    reference_cif = reference_cif or (config.INTERIM_DIR / "pdb" / "6EQE.cif")
    progress_path = progress_path or (config.INTERIM_DIR / "fold_progress.jsonl")

    if max_length:
        n_def = defer_long_candidates(max_length)
        skipped = [(c, s) for c, s in candidates if len(s) > max_length]
        candidates = [(c, s) for c, s in candidates if len(s) <= max_length]
        if skipped:
            print(f"deferred {len(skipped)} candidates over {max_length} aa "
                  f"(longest {max(len(s) for _, s in skipped)} aa); {n_def} rows marked",
                  flush=True)

    todo = [(cid, seq) for cid, seq in candidates if not already_done(cid, out_dir)]
    if limit:
        todo = todo[:limit]
    print(f"{len(candidates)} candidates, {len(candidates)-len(todo)} already folded, "
          f"{len(todo)} to do", flush=True)
    if not todo:
        return []

    t_load = time.monotonic()
    tok, model = load_model()
    print(f"model loaded in {time.monotonic()-t_load:.0f}s", flush=True)

    results: List[FoldResult] = []
    for i, (cid, seq) in enumerate(todo, start=1):
        t0 = time.monotonic()
        try:
            pdb, plddt = fold_one(tok, model, seq)
        except Exception as exc:
            print(f"  [{i}/{len(todo)}] {cid} FAILED: {type(exc).__name__}: {exc}", flush=True)
            continue

        with _TMALIGN_LOCK:
            cif, rmsd, frac = superpose_onto_reference(pdb, reference_cif)

        path = out_dir / f"{cid}.cif"
        path.write_text(cif if cif else "")
        if not cif:
            # Superposition failed: keep the raw prediction rather than losing the fold.
            path = out_dir / f"{cid}.pdb"
            path.write_text(pdb)
        else:
            # The viewer eats PDB, not mmCIF: gemmi's mmCIF omits _entity_poly_seq, so
            # viewers cannot classify the chain as a polymer and refuse to draw a cartoon.
            # The PDB carries the superposed coordinates plus computed HELIX/SHEET records.
            write_viewer_pdb(cif, out_dir / f"{cid}.pdb")

        site = geometry.measure(path)
        el = time.monotonic() - t0
        _persist(cid, path, plddt, frac, rmsd, site)

        with progress_path.open("a") as fh:
            fh.write(json.dumps({"candidate_id": cid, "seconds": round(el, 1),
                                 "plddt": round(plddt, 1), "at": now()}) + "\n")

        results.append(FoldResult(cid, str(path), plddt, frac, rmsd, el))
        print(f"  [{i}/{len(todo)}] {cid} {len(seq)}aa  {el:.0f}s  pLDDT {plddt:.1f}  "
              f"triad {'yes' if site.triad_is_connected else 'NO'}  "
              f"cleft {site.cleft_width_A}", flush=True)

    return results


def write_viewer_pdb(cif_text: str, dest: Path) -> Optional[Path]:
    """Superposed coordinates as PDB, with computed HELIX/SHEET records prepended."""
    import gemmi
    try:
        st = gemmi.read_structure_string(cif_text, format=gemmi.CoorFormat.Mmcif)
        st.setup_entities()
        body = st.make_pdb_string()
    except Exception:
        return None
    ss = secondary_structure(body)
    dest.write_text(_ss_records(ss) + body)
    return dest


def _persist(cid: str, path: Path, plddt: float, frac: Optional[float],
             rmsd: Optional[float], site: geometry.ActiveSite) -> None:
    def _do() -> None:
        with connect() as conn:
            conn.execute(
                "INSERT INTO structures (candidate_id, structure_method, mmcif_path, "
                " plddt_mean, tm_score_to_ispetase, rmsd_ca_to_ispetase_A, "
                " superposition_reference, predicted_at, model_version) "
                "VALUES (?,?,?,?,?,?,?,?,?) "
                "ON CONFLICT(candidate_id) DO UPDATE SET mmcif_path=excluded.mmcif_path, "
                " plddt_mean=excluded.plddt_mean, "
                " tm_score_to_ispetase=excluded.tm_score_to_ispetase, "
                " rmsd_ca_to_ispetase_A=excluded.rmsd_ca_to_ispetase_A, "
                " predicted_at=excluded.predicted_at",
                (cid, "esmfold", path.name, plddt, frac, rmsd,
                 config.ISPETASE_REFERENCE_PDB, now(), config.ESMFOLD_MODEL),
            )
            conn.execute(
                # Every measured column is refreshed on conflict, not just two of them.
                # It used to update only cleft_width_A and the aromatic clamp, so
                # re-running after a fix to any other measurement silently kept the old
                # value -- which is exactly what happened when the oxyanion detection was
                # corrected: the code was right and the stored numbers stayed wrong.
                "INSERT INTO geometry (candidate_id, triad_ser_resnum, triad_asp_resnum, "
                " triad_his_resnum, ser_og_his_ne2_dist_A, his_nd1_asp_od_dist_A, "
                " ser_his_asp_angle_deg, oxyanion_n1_dist_A, oxyanion_n2_dist_A, "
                " oxyanion_n1_resnum, oxyanion_n2_resnum, oxyanion_n2_angle_deg, "
                " cleft_width_A, cleft_depth_A, aromatic_clamp_residues_json) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?) "
                "ON CONFLICT(candidate_id) DO UPDATE SET "
                " triad_ser_resnum=excluded.triad_ser_resnum, "
                " triad_asp_resnum=excluded.triad_asp_resnum, "
                " triad_his_resnum=excluded.triad_his_resnum, "
                " ser_og_his_ne2_dist_A=excluded.ser_og_his_ne2_dist_A, "
                " his_nd1_asp_od_dist_A=excluded.his_nd1_asp_od_dist_A, "
                " ser_his_asp_angle_deg=excluded.ser_his_asp_angle_deg, "
                " oxyanion_n1_dist_A=excluded.oxyanion_n1_dist_A, "
                " oxyanion_n2_dist_A=excluded.oxyanion_n2_dist_A, "
                " oxyanion_n1_resnum=excluded.oxyanion_n1_resnum, "
                " oxyanion_n2_resnum=excluded.oxyanion_n2_resnum, "
                " oxyanion_n2_angle_deg=excluded.oxyanion_n2_angle_deg, "
                " cleft_width_A=excluded.cleft_width_A, "
                " cleft_depth_A=excluded.cleft_depth_A, "
                " aromatic_clamp_residues_json=excluded.aromatic_clamp_residues_json",
                (cid, site.ser_resnum, site.asp_resnum, site.his_resnum,
                 site.ser_og_his_ne2_A, site.his_nd1_asp_od_A, site.ser_his_asp_angle_deg,
                 site.oxyanion_n1_A, site.oxyanion_n2_A,
                 site.oxyanion_n1_resnum, site.oxyanion_n2_resnum,
                 site.oxyanion_n2_angle_deg, site.cleft_width_A,
                 site.cleft_depth_A, json.dumps(site.aromatic_clamp)),
            )
    retry_write(_do)
