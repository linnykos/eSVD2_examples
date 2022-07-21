#!/bin/bash
#$ -N regev_epi_entprog_ni_esvd
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=50G

Rscript --no-save regev_epi_entprog-noninflamed_esvd.R
