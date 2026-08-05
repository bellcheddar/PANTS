"""Does active-site geometry track PET activity, or only fold membership?

The sequence head answered this badly: it separates polyesterases from other folds at
AUC 0.975, and PET-active from PET-inactive polyesterases at chance. Geometry is the one
signal in this project measured off coordinates rather than inherited from somebody's
annotation, so it is the obvious thing to test next, and it is worth testing on exactly
the contrast the embeddings failed.

Structures come from AlphaFold rather than ESMFold: these are characterised enzymes with
UniProt accessions, so a model already exists and folding them again would take hours to
reproduce it. pLDDT is carried through, because a geometric difference that tracks model
confidence is a modelling artefact rather than a finding.
"""

import concurrent.futures as cf
import json
import pathlib
import sys
import urllib.request

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))

from pipeline import config, seqtools
from pipeline.db import connect
from pipeline.recall import pazy
from pipeline.structure import geometry

CACHE = config.INTERIM_DIR / "alphafold"
FEATURES = ["ser_og_his_ne2_A", "his_nd1_asp_od_A", "oxyanion_n1_A", "oxyanion_n2_A",
            "oxyanion_n2_angle_deg", "cleft_width_A", "cleft_depth_A",
            "n_cleft_residues", "n_aromatic_lining"]


API = "https://alphafold.ebi.ac.uk/api/prediction/{acc}"


def fetch(acc: str, dest: pathlib.Path):
    """Resolve the current model URL from the API, then download it."""
    if dest.exists() and dest.stat().st_size > 0:
        return dest
    try:
        with urllib.request.urlopen(API.format(acc=acc), timeout=60) as r:
            meta = json.load(r)
        url = (meta or [{}])[0].get("pdbUrl")
        if not url:
            return None
        urllib.request.urlretrieve(url, dest)
        return dest
    except Exception:
        return None


def mean_plddt(path: pathlib.Path):
    """AlphaFold writes per-residue pLDDT into the B-factor column."""
    vals = [float(l[60:66]) for l in path.read_text().splitlines()
            if l.startswith("ATOM") and l[12:16].strip() == "CA"]
    return sum(vals) / len(vals) if vals else None


if __name__ == "__main__":
    CACHE.mkdir(parents=True, exist_ok=True)
    prots = {f"PAZy:{p['id']}": p for p in pazy.fetch_all()}

    with connect() as c:
        rows = list(c.execute(
            "SELECT enzyme_id, sequence, source_ref FROM characterised_enzymes "
            "WHERE source_ref IN ('PAZy-measured','PAZy-nonPET') AND sequence IS NOT NULL"))

    targets = []
    for eid, seq, src in rows:
        key = eid[:-len("-nonPET")] if src == "PAZy-nonPET" else eid
        p = prots.get(key)
        url = (p or {}).get("alphafold_url")
        if not url:
            continue
        # PAZy stores the EBI *entry page*, not the model. Resolve the coordinate URL
        # from AlphaFold's API rather than composing it: the filename carries a model
        # version, and building "...-model_v4.pdb" by hand returned 404 for every entry
        # once the database moved to v6. A hardcoded version in an external URL is a
        # silent breakage waiting for someone else's release schedule.
        acc = url.rstrip("/").split("/")[-1]
        targets.append((eid, seq, src, acc))
    print(f"{len(targets)} enzymes with an AlphaFold model to fetch", flush=True)

    with cf.ThreadPoolExecutor(max_workers=8) as ex:
        paths = list(ex.map(lambda t: fetch(t[3], CACHE / f"{t[0].replace(':','_')}.pdb"),
                            targets))
    ok = [(t, p) for t, p in zip(targets, paths) if p is not None]
    print(f"downloaded {len(ok)}", flush=True)

    out = []
    for (eid, seq, src, _u), path in ok:
        try:
            site = geometry.measure(path)
        except Exception:
            continue
        if site.ser_resnum is None:
            continue                      # no spatially connected triad found
        rec = {"enzyme_id": eid, "label": 1 if src == "PAZy-measured" else 0,
               "plddt": mean_plddt(path), "seq": seq,
               "triad_connected": site.triad_is_connected}
        for f in FEATURES:
            rec[f] = getattr(site, f, None)
        out.append(rec)

    print(f"measured {len(out)}  "
          f"(PET-active {sum(r['label'] for r in out)}, "
          f"within-family negatives {sum(1 for r in out if not r['label'])})")

    dest = config.ROOT_DIR / "release" / "geometry_vs_activity.json"
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_text(json.dumps(out, indent=2))
    print(f"wrote {dest}")
