#!/bin/bash
#$ -N regev_epi_ta2_ni_esvd
#$ -j y
#$ -o ../../../out/Writeup13/qsub/
#$ -l m_mem_free=50G

Rscript --no-save regev_epi_ta2-noninflamed_esvd.R
