#!/bin/bash
#$ -N cortical_zinbwave
#$ -j y
#$ -o ../../../../out/writeup6b/qsub/
#$ -l m_mem_free=50G
#$ -pe openmpi 4

Rscript --no-save writeup6b_dropseq_humancortical_zinbwave.R
