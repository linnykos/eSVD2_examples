#!/bin/bash
#$ -N sns_layer4
#$ -j y
#$ -o ../../../../out/Writeup13/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_layer4_esvd.R
