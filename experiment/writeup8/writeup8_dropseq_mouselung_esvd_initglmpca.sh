#!/bin/bash
#$ -N lung_esvd_initglmpca
#$ -j y
#$ -o ../../../../out/writeup8/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8_dropseq_mouselung_esvd_initglmpca.R
