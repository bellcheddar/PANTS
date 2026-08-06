#!/usr/bin/env bash
# Push PANTS to the droplet: code, the precomputed database, and the static structures.
#
#   ./deploy.sh              # normal push
#   ./deploy.sh --no-data    # code only, leave the database and structures alone
#
# ALL VALIDATION HAPPENS BEFORE ANYTHING IS COPIED. That ordering is the whole point.
# The chem_sage and chatPDB deploy scripts checked for a missing data/ directory only
# after the code rsync, so a failed deploy left NEW code on disk with the service still
# running the OLD build, and nothing said so. Here, a failure before the first rsync
# leaves the droplet exactly as it was.
set -euo pipefail

DROPLET="${DROPLET_SSH:-root@45.55.102.228}"
APP=pants
ROOT="/opt/$APP"
DOMAIN=pants.mdeller.com
WITH_DATA=1
[ "${1:-}" = "--no-data" ] && WITH_DATA=0

fail() { echo "ABORT: $*" >&2; exit 1; }

# ----------------------------------------------------------------------------------
# Preflight, entirely local, before a single byte moves.
# ----------------------------------------------------------------------------------
echo "==> preflight"
[ -f wsgi.py ] || fail "wsgi.py missing"
[ -f deploy/gunicorn.conf.py ] || fail "deploy/gunicorn.conf.py missing"
[ -f requirements-web.txt ] || fail "requirements-web.txt missing"

for f in release/evaluation_protocol.json release/structure_source_confound.json; do
    [ -f "$f" ] || fail "$f missing -- the Methods page reads its figures from it and \
renders an empty section without it. Run scripts/run_evaluation_protocol.py."
done

if [ "$WITH_DATA" = 1 ]; then
    [ -f pants.db ] || fail "pants.db missing (use --no-data to deploy code only)"
    [ -d app/static/structures ] || fail "app/static/structures missing"
    n_struct=$(find app/static/structures -name '*.pdb' | wc -l | tr -d ' ')
    [ "$n_struct" -gt 0 ] || fail "app/static/structures contains no .pdb files"
    echo "    database $(du -h pants.db | cut -f1), $n_struct structure files"
fi

# The thin venv must be able to import the app. This catches the failure that matters
# most: something in the web path acquiring a dependency the droplet will not have.
if [ -x .venv-web/bin/python ]; then
    .venv-web/bin/python -c "
import wsgi
assert wsgi.app is not None
" >/dev/null 2>&1 || fail "the THIN venv cannot import wsgi:app -- the web layer has \
picked up a dependency the droplet does not have (torch? numpy? biotite?). Fix that \
rather than adding it to requirements-web.txt."
    echo "    thin venv imports wsgi:app cleanly"
else
    echo "    WARNING: no .venv-web locally, skipping the thin-import check"
fi

grep -q "^torch\|^transformers\|^biotite" requirements-web.txt \
    && fail "requirements-web.txt has grown a heavy dependency"

echo "==> pushing code"
rsync -az --delete \
    --exclude '__pycache__' --exclude '*.pyc' --exclude '.DS_Store' \
    app/ "$DROPLET:$ROOT/app/" --exclude 'static/structures'
rsync -az --exclude '__pycache__' --exclude '*.pyc' pipeline/ "$DROPLET:$ROOT/pipeline/"
rsync -az deploy/ "$DROPLET:$ROOT/deploy/"
rsync -az wsgi.py requirements-web.txt "$DROPLET:$ROOT/"

# The small analysis artefacts the pages READ FROM rather than restate. A few kB, and
# they ship with the code because they are what the prose used to be: the Methods page
# reads its AUCs from evaluation_protocol.json, and without the file the paragraph
# renders as nothing at all -- worse than the stale numbers it replaced.
rsync -az --include '*.json' --exclude '*' release/ "$DROPLET:$ROOT/release/"

if [ "$WITH_DATA" = 1 ]; then
    echo "==> pushing data (database + structures)"
    rsync -az pants.db "$DROPLET:$ROOT/pants.db"
    rsync -az --delete --exclude '.DS_Store' app/static/structures/ "$DROPLET:$ROOT/app/static/structures/"
fi

echo "==> installing dependencies, fixing ownership, restarting"
ssh "$DROPLET" bash -euo pipefail <<REMOTE
$ROOT/.venv/bin/pip install -q --upgrade pip
$ROOT/.venv/bin/pip install -q -r $ROOT/requirements-web.txt

# chown EVERY directory, not just the top level. A previous deploy on this box left
# rsync-created subdirectories owned by root, and the service failed to read them.
chown -R $APP:$APP $ROOT

systemctl restart pants-web.service
sleep 2
systemctl is-active --quiet pants-web.service || {
    journalctl -u pants-web.service -n 30 --no-pager
    echo "SERVICE FAILED TO START" >&2
    exit 1
}
REMOTE

# ----------------------------------------------------------------------------------
# http2: certbot rewrites the vhost as `listen 443 ssl;` and drops http2 every time it
# renews or is re-run. Every vhost on this droplet needs it, so it is re-applied on
# each deploy rather than assumed to have survived.
# ----------------------------------------------------------------------------------
echo "==> ensuring http2 on the TLS listeners"
ssh "$DROPLET" bash -euo pipefail <<REMOTE
f=/etc/nginx/sites-available/$APP
if grep -q "listen 443 ssl;" \$f; then
    sed -i 's/listen 443 ssl;/listen 443 ssl http2;/; s/listen \[::\]:443 ssl;/listen [::]:443 ssl http2;/' \$f
    nginx -t && systemctl reload nginx
    echo "    patched"
else
    echo "    already http2 (or TLS not yet configured)"
fi
REMOTE

echo "==> verifying with a real GET"
# A real GET, never curl -I. HEAD can be answered without the response headers Flask
# attaches, so a HEAD check can report this healthy when it is not.
code=$(curl -s -o /dev/null -w '%{http_code}' "https://$DOMAIN/" || echo 000)
if [ "$code" = "200" ]; then
    cc=$(curl -s -D - -o /dev/null "https://$DOMAIN/" | grep -i '^cache-control' || true)
    echo "    https://$DOMAIN/ -> 200"
    echo "    ${cc:-no Cache-Control header on HTML (expected no-cache)}"
else
    echo "    https://$DOMAIN/ -> $code (TLS may not be configured yet; run certbot)"
fi
echo "==> done"
