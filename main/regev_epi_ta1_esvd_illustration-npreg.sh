#!/bin/bash
#$ -N regev_epita1_npreg
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save regev_epi_ta1_esvd_illustration-npreg.R
