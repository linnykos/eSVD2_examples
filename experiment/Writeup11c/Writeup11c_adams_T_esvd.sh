#!/bin/bash
#$ -N adams_T
#$ -j y
#$ -o ../../../../out/Writeup11c/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11c_adams_T_esvd.R
