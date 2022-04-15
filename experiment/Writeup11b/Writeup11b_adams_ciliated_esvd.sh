#!/bin/bash
#$ -N adams_ciliated_esvd
#$ -j y
#$ -o ../../../../out/Writeup11b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup11b_adams_ciliated_esvd.R
