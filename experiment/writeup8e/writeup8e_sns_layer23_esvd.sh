#!/bin/bash
#$ -N sns_layer23_esvd
#$ -j y
#$ -o ../../../../out/writeup8e/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8e_sns_layer23_esvd.R
