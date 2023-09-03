#!/bin/bash
#$ -N regev_epi_umap
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=15G

Rscript --no-save regev_epi-celltype_umap.R
