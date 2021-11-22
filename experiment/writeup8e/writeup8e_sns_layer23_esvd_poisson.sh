#!/bin/bash
#$ -N sns_layer23_esvd_poisson
#$ -j y
#$ -o ../../../../out/writeup8e/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8e_sns_layer23_esvd_poisson.R
