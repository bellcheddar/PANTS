#!/usr/bin/env bash
# One-time provisioning for PANTS on the droplet. Idempotent, but only needed once;
# routine code and data pushes go through deploy.sh.
#
#   ./deploy/provision.sh
#
# Afterwards, run certbot ONCE to obtain the certificate and rewrite the vhost for TLS:
#   ssh root@45.55.102.228 'certbot --nginx -d pants.mdeller.com'
# then re-run this script's http2 patch step (or ./deploy.sh, which repeats it), because
# certbot writes `listen 443 ssl;` without http2 and every vhost on this box needs it.
set -euo pipefail

DROPLET="${DROPLET_SSH:-root@45.55.102.228}"
APP=pants
ROOT="/opt/$APP"
DOMAIN=pants.mdeller.com

echo "==> provisioning $APP on $DROPLET"

ssh "$DROPLET" bash -euo pipefail <<REMOTE
# Service account: no login shell, owns only its own directory.
id -u $APP >/dev/null 2>&1 || useradd --system --home $ROOT --shell /usr/sbin/nologin $APP
mkdir -p $ROOT
chown $APP:$APP $ROOT

# The thin venv. python3-venv is already present (five other apps use it).
if [ ! -x $ROOT/.venv/bin/python ]; then
    python3 -m venv $ROOT/.venv
fi
REMOTE

echo "==> installing systemd unit and nginx site"
scp -q deploy/pants-web.service "$DROPLET:/etc/systemd/system/pants-web.service"
scp -q deploy/nginx.conf "$DROPLET:/etc/nginx/sites-available/$APP"

ssh "$DROPLET" bash -euo pipefail <<REMOTE
ln -sf /etc/nginx/sites-available/$APP /etc/nginx/sites-enabled/$APP
systemctl daemon-reload
systemctl enable pants-web.service
nginx -t && systemctl reload nginx
REMOTE

echo "==> done. Next:"
echo "    ./deploy.sh                                        # push code + data"
echo "    ssh $DROPLET 'certbot --nginx -d $DOMAIN'          # TLS, once"
echo "    ./deploy.sh                                        # re-applies the http2 patch"
