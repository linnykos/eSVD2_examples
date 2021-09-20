#!/bin/bash
#$ -N brain_esvd
#$ -j y
#$ -o ../../../../out/writeup8/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8b_smartseq_mousebrain_esvd.R
