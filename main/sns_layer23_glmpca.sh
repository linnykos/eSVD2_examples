#!/bin/bash
#$ -N sns_layer23_glmpca
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_layer23_glmpca.R
