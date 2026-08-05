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

    @app.after_request
    def _no_heuristic_caching_on_html(response):
        """Stop browsers heuristically caching the HTML.

        Flask sends no Cache-Control on a rendered template. With no explicit header a
        browser is free to apply *heuristic* caching, and it does: it pins the page, and
        with it the `?v=<mtime>` asset URLs embedded in that page. The CSS and JS are then
        cached hard and correctly at those URLs, so a deploy changes the files on disk,
        changes the mtimes, and the browser never learns because it is still reading the
        old HTML. The symptom is a deploy that appears to do nothing.

        Only HTML is marked no-cache. Static assets keep the immutable long-cache header
        from nginx, which is safe precisely because their URLs are content-versioned.

        Verify with a real GET. `curl -I` sends HEAD, which several proxies answer without
        the response headers Flask would attach, so a HEAD check can show this working
        when it is not.
        """
        if response.mimetype == "text/html":
            response.headers["Cache-Control"] = "no-cache, must-revalidate"
        return response

    return app
