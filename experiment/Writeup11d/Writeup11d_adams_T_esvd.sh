#!/bin/bash
#$ -N adams_T
#$ -j y
#$ -o ../../../../out/Writeup11d/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11d_adams_T_esvd.R
