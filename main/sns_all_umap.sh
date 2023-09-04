#!/bin/bash
#$ -N sns_umap
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_all_umap.R
