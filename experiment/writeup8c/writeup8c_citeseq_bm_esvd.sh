#!/bin/bash
#$ -N bm_esvd
#$ -j y
#$ -o ../../../../out/writeup8c/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8c_citeseq_bm_esvd.R
