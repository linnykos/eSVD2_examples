#!/bin/bash
#$ -N regev_epi_ta1_ni_npreg
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=15G

Rscript --no-save regev_epi_ta1-noninflamed_esvd_illustration-npreg.R
