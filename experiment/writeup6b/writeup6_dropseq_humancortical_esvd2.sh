#!/bin/bash
#$ -N cortical_esvd2
#$ -j y
#$ -o ../../../../out/writeup6b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup6_dropseq_humancortical_esvd2.R
