#!/bin/bash
#$ -N adams_T_esvd3
#$ -j y
#$ -o ../../../../out/Writeup12/qsub/
#$ -l m_mem_free=10G

Rscript --no-save adams_T_esvd3.R
