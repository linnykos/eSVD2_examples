#!/bin/bash
#$ -N sns_pseudoreal_extended
#$ -j y
#$ -o ../../../../out/writeup8d/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8d_sns_pseudoreal_esvd_extended.R
