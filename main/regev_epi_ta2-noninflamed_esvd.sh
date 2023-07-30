#!/bin/bash
#$ -N regev_epi_ta2_ni_esvd
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=10G

Rscript --no-save regev_epi_ta2-noninflamed_esvd.R
