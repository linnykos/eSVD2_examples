#!/bin/bash
#$ -N regev_esvd_ta1_inflamed
#$ -j y
#$ -o ../../../../out/Writeup11h/qsub/
#$ -l m_mem_free=75G

Rscript --no-save Writeup11h_regev_esvd_ta1_inflamed.R
