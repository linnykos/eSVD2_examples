#!/bin/bash
#$ -N sns_endothelial
#$ -j y
#$ -o ../../../../out/Writeup13/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_endothelial_esvd.R
