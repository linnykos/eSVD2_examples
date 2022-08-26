#!/bin/bash
#$ -N sns_layer23_esvd_downsampled05
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=50G

Rscript --no-save sns_layer23_esvd_downsampled-0.5.R
