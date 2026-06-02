#!/usr/bin/env bash
# Submit build_database.py as one job per configuration, in the background.
#
# Artifact placement:
#   - Per-config HDF5 outputs go directly into the configuration folder.
#   - Per-run-root ledger and log:
#       <aurora_root>/run_YYYY_MM_DD/build_database_ledger.txt
#       <aurora_root>/run_YYYY_MM_DD/build_database.log
#   - Dispatcher PID + log + queue snapshot:
#       <DISPATCH_DIR>/build_database_dispatch_<ts>.{pid,log,queue.txt}
#     (default DISPATCH_DIR=<aurora_root>/build_database_dispatch)
#   - Nothing is written to the script directory (which is in git).
#
# Usage:
#   ./submit_build_database.sh                     # build + upload, 1 job at a time
#   WORKERS=8 ./submit_build_database.sh           # parallelize across 8 jobs
#   UPLOAD=0 ./submit_build_database.sh            # build only, no upload
#   DRY_RUN=1 ./submit_build_database.sh           # dry-run uploads
#   INCLUDE_RUNS=run_2025_05_25,run_2025_06_01 ./submit_build_database.sh
#   LIMIT=5 ./submit_build_database.sh
#   RETRY_FAILED=1 ./submit_build_database.sh
#
# Tunables:
#   PYTHON         python interpreter (default: python3)
#   AURORA_ROOT    aurora_backup root
#   DISPATCH_DIR   where dispatcher artifacts go (default: $AURORA_ROOT/build_database_dispatch)
#   S3_BUCKET      s3://phy240060/QMCHAMM
#   OSN_ENDPOINT   https://uri.osn.mghpcc.org
#
# Stop a running dispatcher cleanly:
#   kill -INT $(cat <DISPATCH_DIR>/build_database_dispatch_<ts>.pid)

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

PYTHON="${PYTHON:-python3}"
AURORA_ROOT="${AURORA_ROOT:-/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup}"
DISPATCH_DIR="${DISPATCH_DIR:-$AURORA_ROOT/build_database_dispatch}"
S3_BUCKET="${S3_BUCKET:-s3://phy240060/QMCHAMM}"
OSN_ENDPOINT="${OSN_ENDPOINT:-https://uri.osn.mghpcc.org}"
UPLOAD="${UPLOAD:-1}"
DRY_RUN="${DRY_RUN:-0}"
LIMIT="${LIMIT:-0}"
WORKERS="${WORKERS:-1}"
INCLUDE_RUNS="${INCLUDE_RUNS:-}"
RETRY_FAILED="${RETRY_FAILED:-0}"

mkdir -p "$DISPATCH_DIR"
timestamp="$(date +%Y%m%d_%H%M%S)"
log_file="$DISPATCH_DIR/build_database_dispatch_${timestamp}.log"
pid_file="$DISPATCH_DIR/build_database_dispatch_${timestamp}.pid"
queue_file="$DISPATCH_DIR/build_database_dispatch_${timestamp}.queue.txt"

# Keep the script directory pristine: no __pycache__, no logs, no HDF5s.
export PYTHONDONTWRITEBYTECODE=1

if [[ "$UPLOAD" == "1" ]]; then
  read -rp "OSN RW Access Key: " AWS_ACCESS_KEY_ID_INPUT
  read -rp "OSN RW Secret Key: " AWS_SECRET_ACCESS_KEY_INPUT
  echo
  export AWS_ACCESS_KEY_ID="$AWS_ACCESS_KEY_ID_INPUT"
  export AWS_SECRET_ACCESS_KEY="$AWS_SECRET_ACCESS_KEY_INPUT"
  export AWS_DEFAULT_REGION="us-east-1"
  export AWS_REQUEST_CHECKSUM_CALCULATION=when_required
  export AWS_RESPONSE_CHECKSUM_VALIDATION=when_required
fi

# Build common args for the per-config Python invocations.
single_args=(
  --aurora-root "$AURORA_ROOT"
  --s3-bucket "$S3_BUCKET"
  --endpoint "$OSN_ENDPOINT"
)
if [[ "$UPLOAD" == "1" ]]; then single_args+=(--upload); fi
if [[ "$DRY_RUN" == "1" ]]; then single_args+=(--dry-run-upload); fi

list_args=(
  --aurora-root "$AURORA_ROOT"
  --include-runs "$INCLUDE_RUNS"
)
if [[ "$LIMIT" != "0" ]]; then list_args+=(--limit "$LIMIT"); fi
if [[ "$RETRY_FAILED" == "1" ]]; then list_args+=(--retry-failed); fi

cat >>"$log_file" <<EOF
# build_database submission $timestamp
# aurora_root  = $AURORA_ROOT
# dispatch_dir = $DISPATCH_DIR
# s3_bucket    = $S3_BUCKET
# endpoint     = $OSN_ENDPOINT
# upload       = $UPLOAD (dry_run=$DRY_RUN)
# workers      = $WORKERS
# limit        = $LIMIT
# include      = ${INCLUDE_RUNS:-<all>}
# retry        = $RETRY_FAILED
EOF

# Snapshot the work queue using the per-run-root ledgers.
"$PYTHON" "$SCRIPT_DIR/build_database.py" --list-only "${list_args[@]}" >"$queue_file"
n_jobs="$(wc -l <"$queue_file" | tr -d ' ')"
echo "Queued $n_jobs configurations (workers=$WORKERS)"
echo "Dispatch log: $log_file"
echo "Queue file:   $queue_file"
echo "Per-run logs: $AURORA_ROOT/run_YYYY_MM_DD/build_database.log"
echo "Per-run ldg:  $AURORA_ROOT/run_YYYY_MM_DD/build_database_ledger.txt"

if [[ "$n_jobs" == "0" ]]; then
  echo "Nothing to do."
  exit 0
fi

# Dispatch: one --single job per config via xargs -P. The whole dispatcher
# runs under nohup so the user can detach.
nohup bash -c '
  set -euo pipefail
  export PYTHONDONTWRITEBYTECODE=1
  xargs -a "'"$queue_file"'" -P "'"$WORKERS"'" -I {} \
    "'"$PYTHON"'" "'"$SCRIPT_DIR/build_database.py"'" --single "{}" '"${single_args[*]@Q}"'
' >>"$log_file" 2>&1 &
pid=$!
echo "$pid" >"$pid_file"

if [[ "$UPLOAD" == "1" ]]; then
  unset AWS_ACCESS_KEY_ID AWS_SECRET_ACCESS_KEY AWS_DEFAULT_REGION \
        AWS_REQUEST_CHECKSUM_CALCULATION AWS_RESPONSE_CHECKSUM_VALIDATION
fi

echo "PID:    $pid (written to $pid_file)"
echo "Tail:   tail -f \"$log_file\"   # dispatcher output (mostly empty if workers are healthy)"
echo "Stop:   kill -INT $pid           # propagates to xargs children"
