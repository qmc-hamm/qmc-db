#!/usr/bin/env bash
set -euo pipefail

# Sync the schema HDF5 database to OSN using the dedup registry, so nothing is
# ever uploaded twice. Shows the upload set first, asks before sending, then
# re-lists the bucket to confirm the result.
#
#   ./submit_database_sync.sh                      # everything still missing
#   ONLY_SOURCE=ricky_legacy ./submit_database_sync.sh
#   PER_FILE=1 ./submit_database_sync.sh           # slow path, retries per file
#
# The bucket is public-read, so the planning and verification listings need no
# credentials; only the transfer itself does.

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="${ROOT:-/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup}"
S3_BUCKET="${S3_BUCKET:-s3://phy240060/QMCHAMM}"
OSN_ENDPOINT="${OSN_ENDPOINT:-https://uri.osn.mghpcc.org}"
PYTHON="${PYTHON:-/scratch/sgoswam3/ricky_qmchamm/venv/bin/python}"
export AWS_CLI_PYTHON="${AWS_CLI_PYTHON:-/sw/apps/anaconda3/2024.10/bin/python3}"

common=(--root "$ROOT" --bucket "$S3_BUCKET" --endpoint "$OSN_ENDPOINT" --unsigned-list)
[[ -n "${ONLY_SOURCE:-}" ]] && common+=(--only-source "$ONLY_SOURCE")
[[ -n "${PER_FILE:-}" ]] && common+=(--per-file)

timestamp="$(date +%Y%m%d_%H%M%S)"
log="${ROOT%/}/osn_db_sync_${timestamp}.log"
sync_py="$SCRIPT_DIR/sync_database_to_osn.py"

echo "=== plan (dry run) ===" | tee "$log"
"$PYTHON" -u "$sync_py" "${common[@]}" --refresh-remote --dry-run 2>&1 | tee -a "$log"

read -rp "Proceed with the real upload? [y/N] " ok
if [[ "${ok,,}" != "y" ]]; then
  echo "Aborted; nothing uploaded." | tee -a "$log"
  exit 0
fi

read -rp "OSN RW Access Key: " AWS_ACCESS_KEY_ID
read -rsp "OSN RW Secret Key: " AWS_SECRET_ACCESS_KEY
echo
export AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY
export AWS_DEFAULT_REGION="us-east-1"
export AWS_REQUEST_CHECKSUM_CALCULATION=when_required
export AWS_RESPONSE_CHECKSUM_VALIDATION=when_required
trap 'unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_DEFAULT_REGION' EXIT

echo "=== upload ===" | tee -a "$log"
# Reuses the listing just cached by the dry run.
"$PYTHON" -u "$sync_py" "${common[@]}" 2>&1 | tee -a "$log"

echo "=== verify (re-list bucket) ===" | tee -a "$log"
"$PYTHON" -u "$sync_py" "${common[@]}" --refresh-remote --dry-run 2>&1 | tee -a "$log"

echo "Log: $log"
echo "Registry: ${ROOT%/}/upload_registry.csv"
