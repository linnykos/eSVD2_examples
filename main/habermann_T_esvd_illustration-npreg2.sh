#!/bin/bash
#$ -N habermann_T_npreg2
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=15G

Rscript --no-save habermann_T_esvd_illustration-npreg2.R
