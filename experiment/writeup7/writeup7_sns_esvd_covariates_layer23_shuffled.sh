#!/bin/bash
#$ -N sns_esvd_layer23_shuffle
#$ -j y
#$ -o ../../../../out/writeup7/qsub/
#$ -l m_mem_free=100G

Rscript --no-save writeup7_sns_esvd_covariates_layer23_shuffled.R
