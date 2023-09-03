#!/bin/bash
#$ -N sns_layer23
#$ -j y
#$ -o ../../../../out/Writeup13b/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_layer23_esvd.R
