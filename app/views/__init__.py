"""One module per tab, registered onto the app here.

v1 ships Home, Catalogue, Candidate and Methods. Compare, Engineer, Submit and Analysis
are v2 (spec section 13) and get their modules when they land.
"""

from __future__ import annotations

from flask import Flask


def register_views(app: Flask) -> None:
    # Placeholder until the tab modules land in Phase 7. Kept so `PANTS.py serve` and the
    # deploy smoke test work against a real app object from Phase 0 onward.
    @app.route("/healthz")
    def healthz():
        return {"status": "ok", "data_version": app.config["DATA_VERSION"]}
