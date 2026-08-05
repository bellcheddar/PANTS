"""Gunicorn configuration for PANTS.

Sized for the droplet, not for the machine this was developed on. The box has 3.9 GB
shared with AlphaFraud, BoltzMaker, chatPDB, ChemSage and FlexAppeal, so the question is
not "how fast can PANTS be" but "what can PANTS take without harming the neighbours".

Two workers, because the app is IO-bound: it reads a 2.5 MB SQLite file and hands 100 MB
of static structure files to nginx, which serves them directly. Nothing here computes.
"""

bind = "127.0.0.1:8005"
workers = 2
threads = 2
worker_class = "gthread"

# A worker that has been handed a request and gone quiet for two minutes is wedged.
timeout = 120
graceful_timeout = 30
keepalive = 5

# Recycle workers periodically. SQLite connections and Flask's Jinja cache are cheap to
# rebuild, and this bounds any slow leak rather than relying on there not being one.
max_requests = 1000
max_requests_jitter = 100

accesslog = "-"
errorlog = "-"
loglevel = "info"
