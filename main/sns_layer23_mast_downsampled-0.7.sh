#!/bin/bash
#$ -N sns_l23_mast_down_07
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=10G

Rscript --no-save sns_layer23_mast_downsampled-0.9.R
