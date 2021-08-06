#!/bin/bash
#$ -N bm_esvd2
#$ -j y
#$ -o ../../../../out/writeup6b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup6b_citeseq_bm_esvd2.R
