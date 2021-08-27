#!/bin/bash
#$ -N sns_esvd_cov_large2
#$ -j y
#$ -o ../../../../out/writeup7/qsub/
#$ -l m_mem_free=400G

Rscript --no-save writeup7_sns_esvd_covariates_large2.R
