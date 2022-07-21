#!/bin/bash
#$ -N regev_epi-celltype_preprocess
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save regev_epi-celltype_preprocess.R
