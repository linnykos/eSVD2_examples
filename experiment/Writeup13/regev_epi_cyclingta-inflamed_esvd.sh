#!/bin/bash
#$ -N regev_epi_cyclingta_i_esvd
#$ -j y
#$ -o ../../../../out/Writeup13/qsub/
#$ -l m_mem_free=15G

Rscript --no-save regev_epi_cyclingta-inflamed_esvd.R
