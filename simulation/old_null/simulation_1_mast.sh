#!/bin/bash
#$ -N simulation1_mast
#$ -j y
#$ -o ../../../out/simulation/qsub/
#$ -l m_mem_free=10G

Rscript --no-save simulation_1_mast.R
