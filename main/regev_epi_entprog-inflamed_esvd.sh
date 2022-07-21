#!/bin/bash
#$ -N regev_epi_entprog_i_esvd
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=50G

Rscript --no-save regev_epi_entprog-inflamed_esvd.R
