#!/bin/bash
#$ -N bm_esvd
#$ -j y
#$ -o ../../../../out/writeup8/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8_citeseq_bm_esvd_nb.R