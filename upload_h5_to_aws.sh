#!/usr/bin/env bash
set -euo pipefail

# Upload constructed HDF5 files to OSN/S3.
#
# Defaults:
#   SOURCE_DIR   = this script directory
#   S3_BUCKET    = s3://phy240060/QMCHAMM
#   OSN_ENDPOINT = https://uri.osn.mghpcc.org
#
# Example:
#   ./upload_h5_to_aws.sh
#   SOURCE_DIR=/path/to/h5 ./upload_h5_to_aws.sh

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="${SOURCE_DIR:-$SCRIPT_DIR}"
S3_BUCKET="${S3_BUCKET:-s3://phy240060/QMCHAMM}"
OSN_ENDPOINT="${OSN_ENDPOINT:-https://uri.osn.mghpcc.org}"
DO_DRYRUN="${DO_DRYRUN:-1}"

read -rp "Source directory for .h5 upload [${SOURCE_DIR}]: " SOURCE_DIR_INPUT
if [[ -n "${SOURCE_DIR_INPUT}" ]]; then
  SOURCE_DIR="${SOURCE_DIR_INPUT}"
fi
if [[ ! -d "${SOURCE_DIR}" ]]; then
  echo "Source directory does not exist: ${SOURCE_DIR}"
  exit 1
fi

read -rp "OSN RW Access Key: " AWS_ACCESS_KEY_ID_INPUT
read -rp "OSN RW Secret Key: " AWS_SECRET_ACCESS_KEY_INPUT
echo

export AWS_ACCESS_KEY_ID="$AWS_ACCESS_KEY_ID_INPUT"
export AWS_SECRET_ACCESS_KEY="$AWS_SECRET_ACCESS_KEY_INPUT"
export AWS_DEFAULT_REGION="us-east-1"
export AWS_REQUEST_CHECKSUM_CALCULATION=when_required
export AWS_RESPONSE_CHECKSUM_VALIDATION=when_required

shopt -s nullglob
h5_files=("$SOURCE_DIR"/*.h5)
if [[ ${#h5_files[@]} -eq 0 ]]; then
  echo "No .h5 files found in $SOURCE_DIR"
  unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_DEFAULT_REGION
  exit 0
fi

target="${S3_BUCKET%/}/"
echo "Found ${#h5_files[@]} .h5 files"
echo "Target: $target"

timestamp="$(date +%Y%m%d_%H%M%S)"
dry_log="${SOURCE_DIR%/}/osn_h5_upload_dryrun_${timestamp}.log"
run_log="${SOURCE_DIR%/}/osn_h5_upload_${timestamp}.log"

if [[ "$DO_DRYRUN" == "1" ]]; then
  : > "$dry_log"
  echo "Running dry-run upload. Log: $dry_log"
  python3 -m awscli --endpoint-url "$OSN_ENDPOINT" s3 cp "$SOURCE_DIR/" "$target" \
    --recursive --exclude "*" --include "*.h5" --dryrun --no-progress >>"$dry_log" 2>&1
fi

: > "$run_log"
echo "Starting real upload. Log: $run_log"
python3 -m awscli --endpoint-url "$OSN_ENDPOINT" s3 cp "$SOURCE_DIR/" "$target" \
  --recursive --exclude "*" --include "*.h5" --no-progress >>"$run_log" 2>&1

echo "Upload complete."
echo "Uploaded files from: $SOURCE_DIR"
echo "Uploaded to:         $target"
echo "Run log:             $run_log"

unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_DEFAULT_REGION
echo "Credentials unset."
