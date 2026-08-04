"""MGnify acquisition: predicted-protein FASTA for assembly analyses.

Two things this module exists to enforce.

**Only assemblies carry proteins.** MGnify's largest plastic-associated studies are 16S
amplicon ("Targeted Locus"), which contain no protein sequences whatsoever. The single
biggest plastisphere study by sample count, MGYS00001767 with 357 samples, is one of
them. Filtering on `experiment-type == 'assembly'` is not an optimisation, it is the
difference between data and nothing.

**Size is checked before download, not after.** Local free disk is ~62 GB and the whole
v1 footprint is budgeted at ~3.2 GB, so `config.MAX_COLLECTION_GB` caps a run and the
size comes from a HEAD request rather than from watching the disk fill up.
"""

from __future__ import annotations

import gzip
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

from .. import config, http

API = "https://www.ebi.ac.uk/metagenomics/api/v1"
HEADERS = {"Accept": "application/json"}

# Labels MGnify uses for amino-acid CDS files. Both forms appear depending on which
# version of the analysis pipeline produced the study.
CDS_LABELS = ("predicted cds (aa)", "predicted cds with annotation")


@dataclass
class AssemblyFile:
    study: str
    analysis: str
    alias: str
    label: str
    url: str
    size_bytes: Optional[int] = None

    @property
    def size_mb(self) -> Optional[float]:
        return round(self.size_bytes / 1e6, 1) if self.size_bytes else None


def _paged(url: str, params: Optional[dict] = None) -> Iterator[dict]:
    """Walk MGnify's JSON:API pagination."""
    params = dict(params or {})
    params.setdefault("page_size", 100)
    while url:
        resp = http.get(url, params=params, headers=HEADERS)
        if resp.status_code != 200:
            return
        body = resp.json()
        for d in body.get("data", []):
            yield d
        url = (body.get("links") or {}).get("next")
        params = None


def assembly_analyses(study: str, limit: Optional[int] = None) -> List[str]:
    """Analysis accessions of experiment-type 'assembly' for a study."""
    out: List[str] = []
    for d in _paged(f"{API}/studies/{study}/analyses"):
        if (d.get("attributes") or {}).get("experiment-type") == "assembly":
            out.append(d["id"])
            if limit and len(out) >= limit:
                break
    return out


def cds_files(study: str, analysis: str) -> List[AssemblyFile]:
    """Amino-acid CDS downloads for one analysis."""
    resp = http.get(f"{API}/analyses/{analysis}/downloads", headers=HEADERS)
    if resp.status_code != 200:
        return []
    files: List[AssemblyFile] = []
    for it in resp.json().get("data", []):
        a = it.get("attributes", {})
        label = ((a.get("description") or {}).get("label") or "").strip()
        alias = a.get("alias") or ""
        if label.lower() not in CDS_LABELS:
            continue
        # Prefer the unannotated CDS file when both exist: identical sequences, smaller.
        url = (it.get("links") or {}).get("self") or ""
        if url:
            files.append(AssemblyFile(study=study, analysis=analysis, alias=alias,
                                      label=label, url=url))
    return files


def probe_size(f: AssemblyFile) -> Optional[int]:
    """Content-Length via a range request, so nothing large is fetched to learn its size."""
    try:
        resp = http.get(f.url, headers={"Range": "bytes=0-0"}, stream=True)
        rng = resp.headers.get("Content-Range")
        resp.close()
        if rng and "/" in rng:
            total = rng.rsplit("/", 1)[1]
            return int(total) if total.isdigit() else None
        cl = resp.headers.get("Content-Length")
        return int(cl) if cl and cl.isdigit() else None
    except Exception:
        return None


def download(f: AssemblyFile, dest_dir: Path, decompress: bool = True) -> Optional[Path]:
    """Fetch one CDS file. Returns the path to the (decompressed) FASTA."""
    dest_dir.mkdir(parents=True, exist_ok=True)
    gz = dest_dir / f.alias
    out = dest_dir / f.alias[:-3] if f.alias.endswith(".gz") else gz

    if out.exists() and out.stat().st_size > 0:
        return out
    if http.download(f.url, gz) is None:
        return None
    if decompress and gz.suffix == ".gz":
        with gzip.open(gz, "rb") as src, out.open("wb") as dst:
            shutil.copyfileobj(src, dst)
        gz.unlink(missing_ok=True)
        return out
    return gz


def collect(studies: List[str], dest_dir: Path, per_study: int = 2,
            max_total_gb: Optional[float] = None) -> Dict[str, object]:
    """Download CDS FASTA across several studies, stopping before the size cap.

    Returns a report. Nothing is downloaded until its size has been probed and the running
    total checked, so the cap is enforced in advance rather than discovered on a full disk.
    """
    max_total_gb = config.MAX_COLLECTION_GB if max_total_gb is None else max_total_gb
    budget = max_total_gb * 1e9

    report: Dict[str, object] = {
        "downloaded": [], "skipped_over_budget": [], "no_cds": [], "total_bytes": 0,
    }
    total = 0

    for study in studies:
        analyses = assembly_analyses(study, limit=per_study)
        if not analyses:
            report["no_cds"].append(f"{study} (no assembly analyses)")
            continue
        got_any = False
        for an in analyses:
            for f in cds_files(study, an)[:1]:      # one CDS file per analysis
                f.size_bytes = probe_size(f)
                if f.size_bytes and total + f.size_bytes > budget:
                    report["skipped_over_budget"].append(
                        f"{study}/{an} {f.alias} ({f.size_mb} MB)")
                    continue
                path = download(f, dest_dir)
                if path is None:
                    continue
                got_any = True
                total += f.size_bytes or path.stat().st_size
                report["downloaded"].append({
                    "study": study, "analysis": an, "alias": f.alias,
                    "path": str(path), "size_mb": f.size_mb,
                    "bytes_on_disk": path.stat().st_size,
                })
        if not got_any:
            report["no_cds"].append(study)

    report["total_bytes"] = total
    return report


def count_fasta(path: Path) -> Tuple[int, int]:
    """(n_sequences, total_residues) without loading the file into memory."""
    n = res = 0
    with path.open() as fh:
        for line in fh:
            if line.startswith(">"):
                n += 1
            else:
                res += len(line.strip())
    return n, res


# --------------------------------------------------------------------------------------
# Environment labels
# --------------------------------------------------------------------------------------
# Taken from each study's MGnify biome rather than inferred from the filename. The first
# recall run guessed from filename prefixes and left 36 candidates labelled 'unknown',
# which would have silently emptied the Catalogue's environment filter and mis-coloured
# the Home hero plot. Environment is a ranking-relevant axis here, not decoration.
STUDY_BIOME = {
    "MGYS00006036": ("compost", "root:Engineered:Solid waste:Composting"),
    "MGYS00006544": ("marine_plastisphere", "root:Environmental:Aquatic:Marine"),
    "MGYS00004882": ("landfill", "root:Engineered:Solid waste:Landfill"),
    "MGYS00005026": ("compost", "root:Engineered:Solid waste:Composting"),
    "MGYS00004904": ("wastewater", "root:Engineered:Wastewater"),
}

# Which ERZ assembly accession came from which study, so a downloaded FASTA can be
# resolved back to its biome.
ASSEMBLY_STUDY = {
    "ERZ10545954": "MGYS00006036",
    "ERZ21854028": "MGYS00006544",
    "ERZ21854033": "MGYS00006544",
    "ERZ782921":   "MGYS00004882",
    "ERZ794970":   "MGYS00005026",
    "ERZ794971":   "MGYS00005026",
    "ERZ795023":   "MGYS00004904",
}


def environment_for(assembly_or_filename: str) -> str:
    """Environment label for an ERZ accession or a FASTA filename starting with one."""
    for erz, study in ASSEMBLY_STUDY.items():
        if assembly_or_filename.startswith(erz):
            return STUDY_BIOME[study][0]
    return "unknown"


def fetch_biome(study: str) -> Optional[str]:
    """Live biome lookup, for studies not in STUDY_BIOME yet."""
    resp = http.get(f"{API}/studies/{study}", headers=HEADERS)
    if resp.status_code != 200:
        return None
    rel = (resp.json().get("data", {}).get("relationships", {})
           .get("biomes", {}).get("data") or [{}])
    return rel[0].get("id") if rel else None
