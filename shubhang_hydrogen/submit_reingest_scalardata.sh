#!/bin/bash
#SBATCH --job-name=reingest_qmc
#SBATCH --partition=secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --time=04:00:00
#SBATCH --output=reingest_scalardata_%j.out

set -euo pipefail
export OPENBLAS_NUM_THREADS=1

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

python3 reingest_scalardata_h5.py --workers "${SLURM_CPUS_PER_TASK:-8}" "$@"

# After a successful batch, rebuild index hashes for changed files:
# python3 ../ricky_legacy/build_database_index.py --hash --workers 32
