#!/bin/bash
#$ -N sns_layer23_mast
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save sns_layer23_mast.R
