#!/bin/bash
#$ -N lung_esvd
#$ -j y
#$ -o ../../../../out/writeup8b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8b_dropseq_mouselung_esvd.R
