#!/bin/bash
#$ -N sns_layer23_deseq2
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=25G

Rscript --no-save sns_layer23_deseq2.R
