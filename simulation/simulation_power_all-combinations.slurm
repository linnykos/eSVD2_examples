#!/bin/bash

#SBATCH --account=biostat
#SBATCH --partition=all-12c128g
#SBATCH --time=48:00:00

#SBATCH --job-name=esvd-generate
#SBATCH --mem-per-cpu=50gb

#SBATCH --error=/home/users/%u/kzlinlab/projects/eSVD2/out/simulation/slurm/slurm_%j.err
#SBATCH --output=/home/users/%u/kzlinlab/projects/eSVD2/out/simulation/slurm/slurm_%j.out

R CMD BATCH --no-save --no-restore simulation_power_all-combinations.R ../../../out/simulation/slurm/$SLURM_JOB_NAME-$SLURM_JOB_ID.Rout
