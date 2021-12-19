#!/bin/bash
#$ -N sns_pseudoreal_esvd
#$ -j y
#$ -o ../../../../out/writeup8e/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8e_sns_pseudoreal_esvd.R
