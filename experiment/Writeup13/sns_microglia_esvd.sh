#!/bin/bash
#$ -N sns_microglia
#$ -j y
#$ -o ../../../../out/Writeup13/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_microglia_esvd.R
