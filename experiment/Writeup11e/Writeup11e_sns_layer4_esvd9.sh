#!/bin/bash
#$ -N sns_layer4_9
#$ -j y
#$ -o ../../../../out/Writeup11e/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11e_sns_layer4_esvd9.R
