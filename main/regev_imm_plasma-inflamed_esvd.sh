#!/bin/bash
#$ -N regev_imm_plasma_i_esvd
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=30G

Rscript --no-save regev_imm_plasma-inflamed_esvd.R
