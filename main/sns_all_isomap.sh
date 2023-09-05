#!/bin/bash
#$ -N sns_isomap
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_all_isomap.R
