#!/bin/bash
#$ -N sns_layer23_esvd2
#$ -j y
#$ -o ../../../../out/writeup8f/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8f_sns_layer23_esvd2.R
