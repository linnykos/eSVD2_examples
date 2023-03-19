#!/bin/bash
#$ -N regev_imm_macro_ni_esvd
#$ -j y
#$ -o ../../../../out/Writeup13/qsub/
#$ -l m_mem_free=15G

Rscript --no-save regev_imm_macro-noninflamed_esvd.R
