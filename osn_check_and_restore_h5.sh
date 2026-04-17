#!/usr/bin/env bash
set -euo pipefail

# Check and restore HDF5 files from OSN/S3 bucket.
#
# Usage:
#   ./osn_check_and_restore_h5.sh list
#   ./osn_check_and_restore_h5.sh check
#   ./osn_check_and_restore_h5.sh restore-missing
#   ./osn_check_and_restore_h5.sh restore-file <filename.h5>
#
# Defaults:
#   SOURCE_DIR   = this script directory
#   S3_BUCKET    = s3://phy240060/QMCHAMM
#   OSN_ENDPOINT = https://uri.osn.mghpcc.org

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SOURCE_DIR="${SOURCE_DIR:-$SCRIPT_DIR}"
S3_BUCKET="${S3_BUCKET:-s3://phy240060/QMCHAMM}"
OSN_ENDPOINT="${OSN_ENDPOINT:-https://uri.osn.mghpcc.org}"
REMOTE_URI="${S3_BUCKET%/}/"

usage() {
  cat <<EOF
Usage:
  $0 list
  $0 check
  $0 restore-missing
  $0 restore-file <filename.h5>
EOF
  exit 1
}

read -rp "OSN RO Access Key: " AWS_ACCESS_KEY_ID_INPUT
read -rp "OSN RO Secret Key: " AWS_SECRET_ACCESS_KEY_INPUT
echo

export AWS_ACCESS_KEY_ID="$AWS_ACCESS_KEY_ID_INPUT"
export AWS_SECRET_ACCESS_KEY="$AWS_SECRET_ACCESS_KEY_INPUT"
export AWS_DEFAULT_REGION="us-east-1"
export AWS_REQUEST_CHECKSUM_CALCULATION=when_required
export AWS_RESPONSE_CHECKSUM_VALIDATION=when_required

list_remote_h5() {
  python3 -m awscli --endpoint-url "$OSN_ENDPOINT" s3 ls "$REMOTE_URI" | \
    awk '{print $4}' | awk '/\.h5$/ {print $0}'
}

list_local_h5() {
  (cd "$SOURCE_DIR" && ls -1 *.h5 2>/dev/null || true)
}

check_missing() {
  mapfile -t remote_files < <(list_remote_h5)
  mapfile -t local_files < <(list_local_h5)

  if [[ ${#remote_files[@]} -eq 0 ]]; then
    echo "No remote .h5 files found under $REMOTE_URI"
    return 0
  fi

  local tmp_local tmp_remote
  tmp_local="$(mktemp)"
  tmp_remote="$(mktemp)"
  printf "%s\n" "${local_files[@]}" | sort -u > "$tmp_local"
  printf "%s\n" "${remote_files[@]}" | sort -u > "$tmp_remote"

  echo "Local directory: $SOURCE_DIR"
  echo "Remote prefix:   $REMOTE_URI"
  echo
  echo "Missing locally (present remotely):"
  comm -13 "$tmp_local" "$tmp_remote" || true
  echo
  echo "Missing remotely (present locally):"
  comm -23 "$tmp_local" "$tmp_remote" || true

  rm -f "$tmp_local" "$tmp_remote"
}

restore_missing() {
  mapfile -t remote_files < <(list_remote_h5)
  mapfile -t local_files < <(list_local_h5)

  local tmp_local tmp_remote
  tmp_local="$(mktemp)"
  tmp_remote="$(mktemp)"
  printf "%s\n" "${local_files[@]}" | sort -u > "$tmp_local"
  printf "%s\n" "${remote_files[@]}" | sort -u > "$tmp_remote"

  mapfile -t missing_local < <(comm -13 "$tmp_local" "$tmp_remote" || true)
  rm -f "$tmp_local" "$tmp_remote"

  if [[ ${#missing_local[@]} -eq 0 ]]; then
    echo "No missing local .h5 files to restore."
    return 0
  fi

  mkdir -p "$SOURCE_DIR"
  for f in "${missing_local[@]}"; do
    echo "Restoring: $f"
    python3 -m awscli --endpoint-url "$OSN_ENDPOINT" s3 cp "${REMOTE_URI}${f}" "${SOURCE_DIR}/${f}"
  done
}

restore_file() {
  local filename="$1"
  mkdir -p "$SOURCE_DIR"
  echo "Restoring single file: $filename"
  python3 -m awscli --endpoint-url "$OSN_ENDPOINT" s3 cp "${REMOTE_URI}${filename}" "${SOURCE_DIR}/${filename}"
}

cmd="${1:-}"
case "$cmd" in
  list)
    list_remote_h5
    ;;
  check)
    check_missing
    ;;
  restore-missing)
    restore_missing
    ;;
  restore-file)
    [[ $# -eq 2 ]] || usage
    restore_file "$2"
    ;;
  *)
    usage
    ;;
esac

unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_DEFAULT_REGION
echo "Credentials unset."
