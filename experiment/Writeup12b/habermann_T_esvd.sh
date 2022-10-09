#!/bin/bash
#$ -N hab_T_esvd
#$ -j y
#$ -o ../../../../out/Writeup12b/qsub/
#$ -l m_mem_free=10G

Rscript --no-save habermann_T_esvd.R
