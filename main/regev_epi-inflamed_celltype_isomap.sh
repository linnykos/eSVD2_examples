#!/bin/bash
#$ -N regev_epi_inf_isomap
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=15G

Rscript --no-save regev_epi-inflamed_celltype_isomap.R
