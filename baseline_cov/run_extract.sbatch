#!/bin/bash
#SBATCH -J n400_basecov_ex
#SBATCH -p mit_quicktest
#SBATCH -c 1
#SBATCH --mem=8G
#SBATCH -t 00:15:00
#SBATCH -o logs/n400_basecov_ex_%j.out
#SBATCH -e logs/n400_basecov_ex_%j.err
# Uncorrected window + baseline means for kept trials. Reads 4 channels per section.

set -euo pipefail
cd "${SLURM_SUBMIT_DIR:?}"
mkdir -p logs
source /etc/profile.d/00-modulepath.sh
module load miniforge/25.11.0-0
export OMP_NUM_THREADS=1
echo "host $(hostname) job ${SLURM_JOB_ID} start $(date -Is)"
python3 extract_baseline.py
echo "end $(date -Is)"
