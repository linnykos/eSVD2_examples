#!/bin/bash
#$ -N regev_epi_cyclingta_i_esvd
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=50G

Rscript --no-save regev_epi_cyclingta-inflamed_esvd.R
