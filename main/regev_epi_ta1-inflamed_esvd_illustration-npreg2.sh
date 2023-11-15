#!/bin/bash
#$ -N regev_epi_ta1_i_npreg2
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=15G

Rscript --no-save regev_epi_ta1-inflamed_esvd_illustration-npreg2.R
