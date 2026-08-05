import pathlib
import sys, time
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent.parent))
from pathlib import Path
from pipeline import config
from pipeline.recall import mgnify
STUDIES = ["MGYS00006069","MGYS00006006","MGYS00006759","MGYS00002425","MGYS00003619"]
print(f"expanding gut sampling: {len(STUDIES)} studies, up to 10 assemblies each", flush=True)
rep = mgnify.collect(STUDIES, config.RAW_DIR/"mgnify", per_study=10)
tot=0
for d in rep["downloaded"]:
    n,_ = mgnify.count_fasta(Path(d["path"])); tot+=n
print(f"  downloaded {len(rep['downloaded'])} new files, {tot:,} proteins, {rep['total_bytes']/1e6:.0f} MB", flush=True)
if rep["skipped_over_budget"]: print("  over budget:", rep["skipped_over_budget"], flush=True)
# resolve every new assembly to its study so nothing lands as 'unknown'
import json
newmap={}
for study in STUDIES:
    for an in mgnify.assembly_analyses(study, limit=12):
        for f in mgnify.cds_files(study, an):
            erz=f.alias.split("_")[0]
            if erz not in mgnify.ASSEMBLY_STUDY: newmap[erz]=study
Path("/tmp/gut_map2.json").write_text(json.dumps(newmap))
print(f"  {len(newmap)} new assembly accessions to register", flush=True)
