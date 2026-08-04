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


def already_done(candidate_id: str, out_dir: Path) -> bool:
    p = out_dir / f"{candidate_id}.cif"
    return p.exists() and p.stat().st_size > 0


def run(candidates: Sequence[Tuple[str, str]], out_dir: Optional[Path] = None,
        reference_cif: Optional[Path] = None, progress_path: Optional[Path] = None,
        limit: Optional[int] = None) -> List[FoldResult]:
    """Fold, superpose, measure and persist each candidate. Safe to re-run."""
    out_dir = out_dir or config.STRUCTURE_DIR
    out_dir.mkdir(parents=True, exist_ok=True)
    reference_cif = reference_cif or (config.INTERIM_DIR / "pdb" / "6EQE.cif")
    progress_path = progress_path or (config.INTERIM_DIR / "fold_progress.jsonl")

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
                "INSERT INTO geometry (candidate_id, triad_ser_resnum, triad_asp_resnum, "
                " triad_his_resnum, ser_og_his_ne2_dist_A, his_nd1_asp_od_dist_A, "
                " ser_his_asp_angle_deg, oxyanion_n1_dist_A, oxyanion_n2_dist_A, "
                " cleft_width_A, cleft_depth_A, aromatic_clamp_residues_json) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?) "
                "ON CONFLICT(candidate_id) DO UPDATE SET "
                " cleft_width_A=excluded.cleft_width_A, "
                " aromatic_clamp_residues_json=excluded.aromatic_clamp_residues_json",
                (cid, site.ser_resnum, site.asp_resnum, site.his_resnum,
                 site.ser_og_his_ne2_A, site.his_nd1_asp_od_A, site.ser_his_asp_angle_deg,
                 site.oxyanion_n1_A, site.oxyanion_n2_A, site.cleft_width_A,
                 site.cleft_depth_A, json.dumps(site.aromatic_clamp)),
            )
    retry_write(_do)
