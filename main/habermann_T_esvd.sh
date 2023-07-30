#!/bin/bash
#$ -N habermann_T_esvd
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=10G

Rscript --no-save habermann_T_esvd.R
