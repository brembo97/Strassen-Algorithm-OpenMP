#!/bin/sh
#SBATCH --job-name=strassen
#SBATCH --partition=nu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=64
#SBATCH --time=0:30:00

cd $SLURM_SUBMIT_DIR

./strassen 100
