"""Tab modules, registered onto the app here.

v1 ships Home, Catalogue, Candidate, Superpose and Methods. Engineer, Submit and Analysis
are v2 (spec section 13).
"""

from __future__ import annotations

from flask import Flask


def register_views(app: Flask) -> None:
    from .main import bp
    app.register_blueprint(bp)
    from .stats import bp as stats_bp
    app.register_blueprint(stats_bp)
