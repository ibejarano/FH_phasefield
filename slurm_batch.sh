#!/bin/bash
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --mem=4G
#SBATCH --ntasks-per-node=1
#SBATCH --time=24:00:00
#SBATCH --partition=debug

# Cargar entorno de Anaconda
eval "$(conda shell.bash hook)"
conda activate fenicsproject
export OMP_NUM_THREADS=1
NAME=job_${SLURM_JOB_ID}

python run_case.py configs/ps_2d.json --case-name $NAME