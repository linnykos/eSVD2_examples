#!/bin/bash
#$ -N cortical_esvd
#$ -j y
#$ -o ../../../../out/writeup8b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8b_dropseq_humancortical_esvd.R
