#!/bin/bash
#$ -N sns_layer23_npreg2
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=20G

Rscript --no-save sns_layer23_esvd_illustration-npreg2.R
