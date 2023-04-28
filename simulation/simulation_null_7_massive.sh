#!/bin/bash
#$ -N simulation_7
#$ -j y
#$ -o ../../../out/simulation/qsub/
#$ -l m_mem_free=50G

Rscript --no-save simulation_null_7_massive.R
