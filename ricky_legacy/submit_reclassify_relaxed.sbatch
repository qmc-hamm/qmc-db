#!/bin/bash
#SBATCH --job-name=ricky_Sk_relax
#SBATCH --output=/scratch/sgoswam3/ricky_qmchamm/logs/ricky_Sk_relax_%j.out
#SBATCH --error=/scratch/sgoswam3/ricky_qmchamm/logs/ricky_Sk_relax_%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=secondary
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=16G

# Re-label DMC with relaxed min-solid / max-clear-liquid S(k) cuts,
# stamp HDF5 attrs, refresh coverage plots. No S(k) recompute.

set -euo pipefail
mkdir -p /scratch/sgoswam3/ricky_qmchamm/logs

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

PY=/scratch/sgoswam3/ricky_qmchamm/venv/bin/python
SCRIPT=/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup/database_work/qmchamm_database/ricky_legacy/classify_dmc_phase_sk.py
COV=/projects/illinois/grants/qmchamm/shared/shubhang/aurora_backup/HP_dataset_ilnur_coverage

echo "host=$(hostname) date=$(date)"

"$PY" -u "$SCRIPT" --reclassify-relaxed --apply --workers "${SLURM_CPUS_PER_TASK:-16}"

echo "=== refreshing coverage tables + plots ==="
cd "$COV"
"$PY" -u update_coverage_with_sg_database.py
"$PY" -u plot_coverage.py

echo "=== summary ==="
cat /scratch/sgoswam3/ricky_qmchamm/phase_classify/classification_summary.json
echo
echo "finished $(date)"
