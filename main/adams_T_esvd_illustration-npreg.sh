#!/bin/bash
#$ -N adams_T_esvd_illustration-npreg
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save adams_T_esvd_illustration-npreg.R
