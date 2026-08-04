"""The five human gut PET hydrolases (HGMP01 to HGMP05).

Source: Zhang et al., "Identification of a PET hydrolytic enzyme from the human gut
microbiome unveils potential plastic biodegradation in human digestive tract",
Int J Biol Macromol 283 (2024) 137732, PMID 39551294. Sequences come from the authors'
SciDB deposit (`Signalp6_11_pep.fasta`), not from anything reconstructed here.

**The numbering differs between the paper and the deposited data, and matching on the
name would have mislabelled all five.** The paper calls them HGMP01 to HGMP05; the
supplementary BLAST calls the same proteins HGMP03, 04, 06, 07 and 08. The mapping below
is pinned by two independent lines of evidence:

  1. Sequence length. The paper's Table 1 gives 275, 341, 323, 282 and 321 aa. The BLAST
     self-hits give 275, 341, 323, 282 and 320. Every length is distinct, so the
     assignment is 1:1 with no ambiguity.
  2. Homologue count. The paper states "homologue search identified a total of 697
     putative HGMP01-like enzymes". Supplementary HGMP03 returned exactly 697 hits.

Two unrelated numbers agreeing is what makes this a determination rather than a guess.

Why these matter more than their number suggests: **HGMP01 has an optimum of 40 degC and
broad tolerance across pH 7.x**. Every other characterised PET hydrolase in this project
sits at 50 to 60 degC and pH 8 to 9, and even IsPETase (40 degC) has a pH optimum of 9.
HGMP01 is the only measured enzyme in the set whose optimum is close to physiological,
which is the entire premise of PANTS.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional

from .. import config, seqtools
from ..db import connect, now, retry_write

SCIDB_FASTA = config.INTERIM_DIR / "papers" / "scidb" / "Signalp6_11_pep.fasta"
PMID = "39551294"
DOI = "10.1016/j.ijbiomac.2024.137732"


@dataclass
class Hgmp:
    name: str                 # the PAPER's name, which is what anyone reading it will use
    gut_genome_id: str        # identifier in the SciDB deposit and the UHGG database
    suppl_name: str           # the deposit's own, different, name for the same protein
    paper_length: int
    has_dlh: bool
    host_taxonomy: str
    n_homologues: int
    topt_c: Optional[float] = None
    ph_note: Optional[str] = None
    notes: str = ""


HGMPS: List[Hgmp] = [
    Hgmp("HGMP01", "GUT_GENOME238302_00589", "HGMP03", 275, True,
         "Clostridia: Acutalibacteraceae", 697, topt_c=40.0,
         ph_note="Broad tolerance across pH 7.x",
         notes="The active one: hydrolyses PET nanoparticles and outperformed all other "
               "candidates. Shares only ~5% identity with IsPETase. Optimum 40 degC at "
               "near-neutral pH, which is the closest to physiological conditions of any "
               "measured PET hydrolase in this project."),
    Hgmp("HGMP02", "GUT_GENOME243637_00613", "HGMP04", 341, False,
         "Bacteroidia: Porphyromonadaceae", 131),
    Hgmp("HGMP03", "GUT_GENOME137663_00143", "HGMP06", 323, False,
         "Bacteroidia: Porphyromonadaceae", 96),
    Hgmp("HGMP04", "GUT_GENOME171691_00743", "HGMP07", 282, True,
         "Alphaproteobacteria: Rhizobiaceae", 1000),
    Hgmp("HGMP05", "GUT_GENOME244370_00064", "HGMP08", 320, True,
         "Bacteroidia: Porphyromonadaceae", 87),
]


def load(fasta: Optional[Path] = None, tolerance: int = 1) -> Dict[str, object]:
    """Write the five HGMPs into characterised_enzymes as measured-activity positives.

    Every sequence is length-checked against the paper's own table before it is written.
    That check is not ceremonial: it caught a transposition between the similar ids
    GUT_GENOME244370_00064 and GUT_GENOME244064_00699 during this curation, which would
    otherwise have stored the wrong protein under HGMP05.
    """
    seqs = seqtools.read_fasta(fasta or SCIDB_FASTA)
    report: Dict[str, object] = {"loaded": [], "problems": []}

    with connect() as conn:
        for h in HGMPS:
            seq = seqs.get(h.gut_genome_id)
            if seq is None:
                report["problems"].append(f"{h.name}: {h.gut_genome_id} absent from FASTA")
                continue
            if abs(len(seq) - h.paper_length) > tolerance:
                report["problems"].append(
                    f"{h.name}: {len(seq)} aa but the paper says {h.paper_length}")
                continue
            if not set(seq) <= set("ACDEFGHIKLMNPQRSTVWY"):
                report["problems"].append(f"{h.name}: non-standard residues")
                continue

            conn.execute(
                "INSERT INTO characterised_enzymes "
                "(enzyme_id, uniprot, organism, family, sequence, seq_length, "
                " is_positive, is_negative, is_near_miss, topt_c, "
                " activity_substrate_notes, source_ref, added_at) "
                "VALUES (?,?,?,?,?,?,1,0,0,?,?,?,?) "
                "ON CONFLICT(enzyme_id) DO UPDATE SET sequence=excluded.sequence, "
                " seq_length=excluded.seq_length, topt_c=excluded.topt_c, "
                " activity_substrate_notes=excluded.activity_substrate_notes",
                (h.name, None, h.host_taxonomy, "petase_like", seq, len(seq), h.topt_c,
                 (f"Human gut PET hydrolase, measured. {h.notes} "
                  f"Deposit id {h.gut_genome_id} (the paper's {h.name} is {h.suppl_name} "
                  f"in the supplementary). {h.n_homologues} homologues in UHGP-100. "
                  f"DLH domain: {'yes' if h.has_dlh else 'no'}. PMID {PMID}.").strip(),
                 "HGMP-measured", now()),
            )
            report["loaded"].append((h.name, h.gut_genome_id, len(seq)))

        # Measured optima, with the citation attached, same as the UniProt extraction.
        for h in HGMPS:
            if h.topt_c is None:
                continue
            conn.execute(
                "DELETE FROM activity_measurements WHERE enzyme_id=? AND extraction_confidence='paper'",
                (h.name,))
            conn.execute(
                "INSERT INTO activity_measurements "
                "(enzyme_id, substrate_form, temperature_c, product_measured, "
                " parameter_type, rate_value, rate_units, raw_text, evidence_code, "
                " comparable_group_id, source_doi, extracted_at, extraction_confidence) "
                "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?)",
                (h.name, "nanoparticle", h.topt_c, "PET monomers", "topt", h.topt_c,
                 "degC",
                 f"Optimum temperature {h.topt_c:.0f} degC on PET nanoparticles. {h.ph_note or ''}",
                 "ECO:0000269", "topt:pet_nanoparticle",
                 f"PMID:{PMID};doi:{DOI}", now(), "paper"),
            )

    def _src() -> None:
        with connect() as conn:
            conn.execute(
                "INSERT INTO data_sources(name, version, retrieved_at, n_records, license, source_url) "
                "VALUES (?,?,?,?,?,?) ON CONFLICT(name) DO UPDATE SET "
                "retrieved_at=excluded.retrieved_at, n_records=excluded.n_records",
                ("HGMP-SciDB", f"PMID {PMID}", now(), len(report["loaded"]),
                 "see publication", f"https://doi.org/{DOI}"))
    retry_write(_src)
    return report
