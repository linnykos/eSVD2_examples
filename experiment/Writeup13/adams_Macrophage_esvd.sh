#!/bin/bash
#$ -N adams_Macrophage_esvd
#$ -j y
#$ -o ../../../../out/Writeup13/qsub/
#$ -l m_mem_free=75G

Rscript --no-save adams_Macrophage_esvd.R
