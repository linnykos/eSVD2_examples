#!/bin/bash
#$ -N sns_layer4_esvd
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=100G

Rscript --no-save sns_layer4_esvd.R
