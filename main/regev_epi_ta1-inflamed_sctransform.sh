#!/bin/bash
#$ -N regev_epi_ta1_i_sctransform
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=150G

Rscript --no-save regev_epi_ta1-inflamed_sctransform.R
