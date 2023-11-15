#!/bin/bash
#$ -N adams_T_esvd_illustration-npreg2
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=20G

Rscript --no-save adams_T_esvd_illustration-npreg2.R
