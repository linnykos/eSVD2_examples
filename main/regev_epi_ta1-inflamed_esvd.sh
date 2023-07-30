#!/bin/bash
#$ -N regev_epi_ta1_i_esvd
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=10G

Rscript --no-save regev_epi_ta1-inflamed_esvd.R
