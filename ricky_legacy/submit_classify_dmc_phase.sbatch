#!/bin/bash
#SBATCH --job-name=ricky_dmc_Sk_phase
#SBATCH --output=/scratch/sgoswam3/ricky_qmchamm/logs/ricky_dmc_Sk_phase_%j.out
#SBATCH --error=/scratch/sgoswam3/ricky_qmchamm/logs/ricky_dmc_Sk_phase_%j.err
#SBATCH --time=04:00:00
#SBATCH --partition=secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=32G

# Full S(k) sweep of all 17,528 Ricky DMC configs vs known solid/liquid refs.
# Writes phase attrs (solid/liquid/ambiguous) and refreshes coverage plots.

set -euo pipefail
mkdir -p /scratch/sgoswam3/ricky_qmchamm/logs /scratch/sgoswam3/ricky_qmchamm/phase_classify

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

PY=/scratch/sgoswam3/ricky_qmchamm/venv/bin/python
SCRIPT=/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup/database_work/qmchamm_database/ricky_legacy/classify_dmc_phase_sk.py
COV=/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup/HP_dataset_ilnur_coverage

echo "host=$(hostname)  date=$(date)  cpus=${SLURM_CPUS_PER_TASK:-16}"
echo "python=$PY"

"$PY" -u "$SCRIPT" --workers "${SLURM_CPUS_PER_TASK:-16}" --apply

echo "=== refreshing coverage tables + plots ==="
cd "$COV"
"$PY" -u update_coverage_with_sg_database.py
"$PY" -u plot_coverage.py

echo "=== summary ==="
cat /scratch/sgoswam3/ricky_qmchamm/phase_classify/classification_summary.json
echo
echo "finished $(date)"
