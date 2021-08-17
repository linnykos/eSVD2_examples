#!/bin/bash
#$ -N sns_esvd
#$ -j y
#$ -o ../../../../out/writeup7/qsub/
#$ -l m_mem_free=400G

Rscript --no-save writeup7_sns_esvd.R
