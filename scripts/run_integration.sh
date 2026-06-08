#!/bin/bash
#SBATCH --job-name=integration
#SBATCH --partition=small
#SBATCH --mem=80G
#SBATCH --cpus-per-task=4
#SBATCH --time=02:00:00
#SBATCH --output=/scratch/svc_td_compbio/users/MaDeBa/figures/logs/integration_%j.log

# Create directory for the log file if it doesn't exist yet
mkdir -p /scratch/svc_td_compbio/users/MaDeBa/figures/logs

# Clean environment
module purge
module load R/4.4.1

# Load and activate conda
source ~/miniforge3/etc/profile.d/conda.sh
conda activate scvi-env-310

# Run Rscript with library overrides
LD_PRELOAD=$CONDA_PREFIX/lib/libstdc++.so.6 \
LD_LIBRARY_PATH=$CONDA_PREFIX/lib:$LD_LIBRARY_PATH \
Rscript /scratch/svc_td_compbio/users/MaDeBa/scripts/analysis_integration_RNA_MDB.R
