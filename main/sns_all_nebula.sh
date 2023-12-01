#!/bin/bash
#$ -N sns_nebula
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=50G

Rscript --no-save sns_all_nebula.R
