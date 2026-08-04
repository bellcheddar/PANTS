"""ESM-2 embeddings: frozen encoder, mean-pooled, CPU.

Spec section 5.5: freeze ESM-2 t12-35M and train shallow heads. Do NOT fine-tune end to
end, because with order 10^2 real positives an end-to-end fine-tune memorises. Everything
here is inference only; nothing in this module has gradients.

Why t12-35M and not 650M: the same weights have to run inside a Flask worker on a 3.8 GB
droplet shared with five other applications (spec section 3.1, PLAN_v1.md section 5).
Embedding offline with 650M and serving with 35M would mean the live `/submit` path scored
sequences in a different representation from the catalogue, which is worse than the
accuracy it would buy.

CPU rather than MPS, deliberately: 35M parameters is small enough that the transfer
overhead and MPS's dtype quirks are not worth it, and the offline run is minutes either way.
"""

from __future__ import annotations

import time
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from .. import config

# Longer than any single-domain alpha/beta hydrolase we keep after filtering (203-450 aa).
# Truncation is reported rather than silent: a truncated embedding is a different object
# from a complete one and should not quietly enter the training set.
MAX_LENGTH = 1024


def load_model(model_name: Optional[str] = None):
    """Tokenizer and model, on CPU, in eval mode."""
    import torch
    from transformers import AutoModel, AutoTokenizer

    name = model_name or config.ESM_MODEL
    tok = AutoTokenizer.from_pretrained(name)
    model = AutoModel.from_pretrained(name)
    model.eval()
    torch.set_grad_enabled(False)
    return tok, model


def embed(records: Sequence[Tuple[str, str]], batch_size: int = 8,
          model_name: Optional[str] = None, progress_every: int = 200
          ) -> Tuple[List[str], np.ndarray, Dict[str, object]]:
    """Mean-pooled per-sequence embeddings.

    Returns (ids, matrix of shape (n, dim), report). Pooling excludes padding AND the
    special CLS/EOS tokens: including padding would make an embedding depend on what else
    happened to be in its batch, which is a genuinely nasty bug because the vectors still
    look reasonable.
    """
    import torch

    tok, model = load_model(model_name)
    ids: List[str] = []
    vecs: List[np.ndarray] = []
    truncated: List[str] = []
    t0 = time.monotonic()

    for start in range(0, len(records), batch_size):
        chunk = records[start:start + batch_size]
        for sid, seq in chunk:
            if len(seq) > MAX_LENGTH - 2:
                truncated.append(sid)

        batch_ids = [c[0] for c in chunk]
        batch_seqs = [c[1][: MAX_LENGTH - 2] for c in chunk]
        enc = tok(batch_seqs, return_tensors="pt", padding=True, truncation=True,
                  max_length=MAX_LENGTH)
        out = model(**enc).last_hidden_state          # (B, L, D)

        mask = enc["attention_mask"].clone()
        # Drop CLS (first real token) and EOS (last real token of each sequence).
        mask[:, 0] = 0
        lengths = enc["attention_mask"].sum(dim=1)
        for i, L in enumerate(lengths.tolist()):
            mask[i, L - 1] = 0

        m = mask.unsqueeze(-1).float()
        pooled = (out * m).sum(dim=1) / m.sum(dim=1).clamp(min=1)
        vecs.append(pooled.numpy().astype(np.float32))
        ids.extend(batch_ids)

        if progress_every and (start // batch_size) % max(1, progress_every // batch_size) == 0:
            done = min(start + batch_size, len(records))
            print(f"  embedded {done}/{len(records)}  ({time.monotonic()-t0:.0f}s)", flush=True)

    matrix = np.vstack(vecs) if vecs else np.zeros((0, config.ESM_EMBED_DIM), dtype=np.float32)
    report = {
        "n": len(ids),
        "dim": int(matrix.shape[1]) if matrix.size else 0,
        "seconds": round(time.monotonic() - t0, 1),
        "seq_per_sec": round(len(ids) / max(1e-9, time.monotonic() - t0), 1),
        "truncated": truncated,
        "model": model_name or config.ESM_MODEL,
    }
    return ids, matrix, report


def save(ids: Sequence[str], matrix: np.ndarray, path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, ids=np.array(ids, dtype=object), embeddings=matrix)
    return path


def load(path: Path) -> Tuple[List[str], np.ndarray]:
    d = np.load(path, allow_pickle=True)
    return list(d["ids"]), d["embeddings"]
