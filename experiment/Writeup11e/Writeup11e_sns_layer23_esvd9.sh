#!/bin/bash
#$ -N sns_layer23_9
#$ -j y
#$ -o ../../../../out/Writeup11e/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11e_sns_layer23_esvd9.R
