#!/bin/bash
#$ -N sns_deseq2_plot
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=20G

Rscript --no-save sns_all_volcano-ggrepel.R
