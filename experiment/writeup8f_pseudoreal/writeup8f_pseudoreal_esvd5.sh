#!/bin/bash
#$ -N pseudoreal_esvd5
#$ -j y
#$ -o ../../../../out/writeup8f/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8f_pseudoreal_esvd5.R