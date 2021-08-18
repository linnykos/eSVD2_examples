#!/bin/bash
#$ -N sns_esvd_cov
#$ -j y
#$ -o ../../../../out/writeup7/qsub/
#$ -l m_mem_free=200G

Rscript --no-save writeup7_sns_esvd_covariates.R
