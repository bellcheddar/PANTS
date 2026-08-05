# Pipeline drivers

Thin scripts that drive the `pipeline/` package. The library holds the logic and the
tests; these hold the sequence of calls that produced the numbers in the project README,
which is why they belong in the repository rather than in a scratch directory.

Run from the repository root with the pipeline virtual environment:

```bash
.venv/bin/python scripts/<name>.py
```

| Script | What it does |
|---|---|
| `curate_ec.py` | Harvests EC 3.1.1.101 from UniProt and extracts measured kinetics |
| `more_gut.py` | Downloads additional MGnify gut assemblies, size-checked first |
| `run_recall.py` | Recall over a set of metagenome FASTA files |
| `run_gut2.py` | Recall over the gut collection; **resumable** via `data/interim/gut_done.txt` |
| `run_embed.py` | ESM-2 embedding of the training set and candidates |
| `run_fold.py` | One-shot ESMFold pass over the current candidate list |
| `fold_drain.py` | **Continuously** folds whatever has no structure yet, and keeps up with recall |
| `run_train.py` | The three-scheme labelling experiment (naive / evidence-only / PU) |
| `retest.py` | Composition baseline and head on the measured-positive set |
| `nearmiss.py` | The within-family contrast: measured positives against near misses |

Two things worth knowing before running any of them.

**The long ones are resumable and should stay that way.** `fold_drain.py` skips
candidates that already have a structure file; `run_gut2.py` skips files listed in
`gut_done.txt`. Both were interrupted several times during the first full run (a length
cap added mid-flight, an external stop, a crash on an empty assembly) and lost nothing but
the work in flight.

**Do not run `fold_drain.py` alongside a recall scan.** ESMFold and MMseqs2 compete for
the same cores: a 287 aa fold took 2,696 s under contention against roughly 110 s with the
machine to itself. `fold_drain.py` waits for `run_fold.py` to exit, but it does not know
about recall.
