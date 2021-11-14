#!/bin/bash
#$ -N sns_layer23_esvd_extended
#$ -j y
#$ -o ../../../../out/writeup8d/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8d_sns_layer23_esvd_extended.R
