import pathlib
import sys
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from pipeline import uniprot
from pipeline.db import connect, now
from pipeline.db.manifest import stage_manifest
from pipeline.harmonise import activity

with stage_manifest("harmonise", label="ec-3.1.1.101") as man:
    ents = list(uniprot.search(f"ec:{activity.EC_PETASE}", max_results=600))
    print(f"EC {activity.EC_PETASE}: {len(ents)} entries ({sum(e.reviewed for e in ents)} reviewed)", flush=True)

    with connect() as c:
        known = {r[0]: r[1] for r in c.execute(
            "SELECT uniprot, enzyme_id FROM characterised_enzymes WHERE uniprot IS NOT NULL")}

    added = upgraded = 0
    with connect() as c:
        for e in ents:
            if not e.is_plausible_protein:
                continue
            eid = known.get(e.accession)
            tier = "EC-3.1.1.101-reviewed" if e.reviewed else "EC-3.1.1.101-unreviewed"
            if eid:
                c.execute("UPDATE characterised_enzymes SET source_ref=?, is_positive=1 WHERE enzyme_id=?", (tier, eid))
                upgraded += 1
            else:
                eid = f"EC:{e.accession}"
                c.execute(
                    "INSERT INTO characterised_enzymes (enzyme_id, uniprot, organism, family, "
                    " sequence, seq_length, is_positive, is_negative, is_near_miss, "
                    " taxonomy_lineage, activity_substrate_notes, source_ref, added_at) "
                    "VALUES (?,?,?,?,?,?,1,0,0,?,?,?,?) ON CONFLICT(enzyme_id) DO NOTHING",
                    (eid, e.accession, e.organism, "petase_like", e.sequence, e.length,
                     e.lineage,
                     f"EC {activity.EC_PETASE} (poly(ethylene terephthalate) hydrolase) assigned by UniProt"
                     + (", Swiss-Prot reviewed" if e.reviewed else ""),
                     tier, now()))
                added += 1
                known[e.accession] = eid
    print(f"positives: {added} added, {upgraded} upgraded to measured-activity tier", flush=True)

    # Kinetics only for the reviewed ones: unreviewed entries carry no curated
    # biophysicochemical block, so fetching all 459 would be ~450 wasted requests.
    reviewed = [e for e in ents if e.reviewed]
    print(f"extracting measurements from {len(reviewed)} reviewed entries...", flush=True)
    all_m = []
    for e in reviewed:
        ms, meta = activity.extract(e.accession, known.get(e.accession, f"EC:{e.accession}"))
        all_m.extend(ms)
    n = activity.write(all_m)
    man.counts(n_input=len(ents), n_output=n, n_discarded=len(ents)-len(reviewed))
    print(f"wrote {n} measurements", flush=True)
