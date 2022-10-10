#!/bin/bash
#$ -N hab_T_esvd_old
#$ -j y
#$ -o ../../../../out/Writeup12b/qsub/
#$ -l m_mem_free=10G

Rscript --no-save habermann_T_esvd_old-v057.R
