#!/bin/bash
#$ -N sns_layer4_esvd3_initialize
#$ -j y
#$ -o ../../../../out/Writeup10/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup10_sns_layer4_esvd3_initialize.R
