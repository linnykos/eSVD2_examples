#!/bin/bash
#$ -N adams_epithelial_esvd
#$ -j y
#$ -o ../../../../out/Writeup11b/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11b_adams_epithelial_esvd.R