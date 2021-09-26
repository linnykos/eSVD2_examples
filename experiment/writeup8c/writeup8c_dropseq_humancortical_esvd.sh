#!/bin/bash
#$ -N cortical_esvd
#$ -j y
#$ -o ../../../../out/writeup8c/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8c_dropseq_humancortical_esvd.R
