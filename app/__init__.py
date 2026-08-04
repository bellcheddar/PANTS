"""Flask application factory.

The web layer reads precomputed SQLite and serves static mmCIF. It must stay light: the
droplet has 3.8 GB shared with five other apps, so torch and transformers are never
imported here in v1 (see requirements-web.txt and PLAN_v1.md section 5).
"""

from __future__ import annotations

from flask import Flask

from pipeline import config


def create_app() -> Flask:
    app = Flask(
        __name__,
        static_folder=str(config.STATIC_DIR),
        template_folder=str(config.APP_DIR / "templates"),
    )
    app.config["DATA_VERSION"] = config.DATA_VERSION
    app.config["SERVER_NAME_DISPLAY"] = config.SERVER_NAME

    from .views import register_views
    register_views(app)

    return app
