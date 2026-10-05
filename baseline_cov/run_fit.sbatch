#!/bin/bash
#SBATCH -J n400_basecov_lmm
#SBATCH -p mit_preemptable
#SBATCH --requeue
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 04:00:00
#SBATCH -o logs/n400_basecov_lmm_%j.out
#SBATCH -e logs/n400_basecov_lmm_%j.err
# Baseline-covariate vs subtraction on the revision core model.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
mkdir -p logs
source /etc/profile.d/00-modulepath.sh
module load R/4.5.2
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
export OPENBLAS_NUM_THREADS="${SLURM_CPUS_PER_TASK:-8}"
echo "host $(hostname) job ${SLURM_JOB_ID} start $(date -Is)"
Rscript fit_basecov.R "${ROI:-cp}"
echo "end $(date -Is)"
