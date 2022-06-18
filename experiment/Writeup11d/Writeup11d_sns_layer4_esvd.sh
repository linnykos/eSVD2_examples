#!/bin/bash
#$ -N sns_layer4
#$ -j y
#$ -o ../../../../out/Writeup11d/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11d_sns_layer4_esvd.R
