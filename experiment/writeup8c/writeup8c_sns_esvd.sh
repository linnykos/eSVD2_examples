#!/bin/bash
#$ -N sns_esvd
#$ -j y
#$ -o ../../../../out/writeup8c/qsub/
#$ -l m_mem_free=150G

Rscript --no-save writeup8c_sns_esvd.R