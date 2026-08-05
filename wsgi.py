"""WSGI entry point for gunicorn.

Deliberately trivial: the app factory is in app/__init__.py and this exists only so the
systemd unit has a stable `wsgi:app` target.
"""

from app import create_app

app = create_app()
