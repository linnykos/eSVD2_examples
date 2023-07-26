#!/bin/bash
#$ -N sns_endothelial_esvd
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_endothelial_esvd.R
