#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=mast21
#SBATCH --mem-per-cpu=50gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/eSVD2/out/simulation/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/eSVD2/out/simulation/slurm/slurm_%j.out

Rscript --no-save --no-restore simulation_power_mast_all-combinations.R 2 1 > ../../../out/simulation/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.Rout
