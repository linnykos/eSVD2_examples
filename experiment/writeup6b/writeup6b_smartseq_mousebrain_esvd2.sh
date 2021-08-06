#!/bin/bash
#$ -N brain_esvd2
#$ -j y
#$ -o ../../../../out/writeup6b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup6b_smartseq_mousebrain_esvd2.R
