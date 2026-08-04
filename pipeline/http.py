"""Shared retrying HTTP session for the external services (UniProt, RCSB, MGnify).

Ported from AlphaFraud's http.py, which has been hitting public bioinformatics APIs in
production since 2026-07-14. One session gives connection reuse; tenacity gives polite
backoff on the transient 5xx and rate-limit responses these APIs occasionally return.
"""

from __future__ import annotations

import os
import threading
from typing import Any, Optional

import requests
from tenacity import (
    retry,
    retry_if_exception_type,
    stop_after_attempt,
    wait_exponential,
)

from . import config

_SESSION = requests.Session()
_SESSION.headers.update({"User-Agent": "PANTS/0.1 (+https://pants.mdeller.com)"})

# Retry only on connection errors and the handful of server-side statuses worth retrying.
_RETRYABLE_STATUS = {429, 500, 502, 503, 504}


class RetryableStatus(Exception):
    """Raised for a retryable HTTP status so tenacity re-attempts the request."""


_retry = retry(
    reraise=True,
    stop=stop_after_attempt(config.HTTP_MAX_RETRIES),
    wait=wait_exponential(multiplier=1, min=1, max=30),
    retry=retry_if_exception_type((requests.ConnectionError, requests.Timeout, RetryableStatus)),
)


@_retry
def get(url: str, **kwargs) -> requests.Response:
    resp = _SESSION.get(url, timeout=config.HTTP_TIMEOUT, **kwargs)
    if resp.status_code in _RETRYABLE_STATUS:
        raise RetryableStatus(f"{resp.status_code} for GET {url}")
    return resp


@_retry
def post_json(url: str, payload: dict[str, Any], **kwargs) -> requests.Response:
    resp = _SESSION.post(url, json=payload, timeout=config.HTTP_TIMEOUT, **kwargs)
    if resp.status_code in _RETRYABLE_STATUS:
        raise RetryableStatus(f"{resp.status_code} for POST {url}")
    return resp


def get_json(url: str, **kwargs) -> Optional[Any]:
    """GET returning parsed JSON, or None on 404."""
    resp = get(url, **kwargs)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    return resp.json()


def graphql(url: str, query: str, variables: Optional[dict] = None) -> dict:
    """POST a GraphQL query; return the `data` object (raises on transport/GraphQL error)."""
    resp = post_json(url, {"query": query, "variables": variables or {}})
    resp.raise_for_status()
    body = resp.json()
    if body.get("errors"):
        raise RuntimeError(f"GraphQL errors: {body['errors']}")
    return body.get("data", {})


def download(url: str, dest, skip_if_exists: bool = True):
    """Stream a URL to `dest` (a Path). Returns dest, or None on 404."""
    from pathlib import Path

    dest = Path(dest)
    if skip_if_exists and dest.exists() and dest.stat().st_size > 0:
        return dest
    resp = get(url, stream=True)
    if resp.status_code == 404:
        return None
    resp.raise_for_status()
    dest.parent.mkdir(parents=True, exist_ok=True)
    # Unique temp name per writer so parallel workers downloading the same file never
    # clobber each other's partial write.
    tmp = dest.with_suffix(f"{dest.suffix}.{os.getpid()}.{threading.get_ident()}.part")
    try:
        with open(tmp, "wb") as fh:
            for chunk in resp.iter_content(chunk_size=1 << 16):
                fh.write(chunk)
        tmp.replace(dest)
    finally:
        if tmp.exists():
            tmp.unlink(missing_ok=True)
    return dest
