#!/bin/bash
#$ -N sns_layer4_esvd_resample
#$ -j y
#$ -o ../../../../out/Writeup11/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11_sns_layer4_esvd_resampling.R
