#!/bin/bash
#$ -N adams_T_esvd
#$ -j y
#$ -o ../../../../out/Writeup11b/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11b_adams_T_esvd.R
