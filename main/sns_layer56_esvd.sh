#!/bin/bash
#$ -N sns_layer56_esvd
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=50G

Rscript --no-save sns_layer56_esvd.R
