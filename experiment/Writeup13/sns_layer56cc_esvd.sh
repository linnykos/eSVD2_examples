#!/bin/bash
#$ -N sns_layer56cc
#$ -j y
#$ -o ../../../../out/Writeup13/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_layer56cc_esvd.R
