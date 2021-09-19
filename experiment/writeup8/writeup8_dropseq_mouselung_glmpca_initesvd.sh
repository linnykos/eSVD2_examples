#!/bin/bash
#$ -N lung_glmpca_initesvd
#$ -j y
#$ -o ../../../../out/writeup8/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8_dropseq_mouselung_glmpca_initesvd.R
