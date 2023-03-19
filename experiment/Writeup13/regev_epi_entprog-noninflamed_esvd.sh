#!/bin/bash
#$ -N regev_epi_entprog_ni_esvd
#$ -j y
#$ -o ../../../../out/Writeup13/qsub/
#$ -l m_mem_free=15G

Rscript --no-save regev_epi_entprog-noninflamed_esvd.R
