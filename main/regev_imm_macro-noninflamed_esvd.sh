#!/bin/bash
#$ -N regev_imm_macro_ni_esvd
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=75G

Rscript --no-save regev_imm_macro-noninflamed_esvd.R
