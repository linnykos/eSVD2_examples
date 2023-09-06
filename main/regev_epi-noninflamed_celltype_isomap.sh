#!/bin/bash
#$ -N regev_epi_noninf_isomap
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=15G

Rscript --no-save regev_epi-noninflamed_celltype_isomap.R
