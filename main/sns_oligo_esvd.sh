#!/bin/bash
#$ -N sns_oligo_esvd
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save sns_oligo_esvd.R
