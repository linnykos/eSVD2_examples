#!/bin/bash
#$ -N hab_Macrophage_esvd
#$ -j y
#$ -o ../../../../out/Writeup13/qsub/
#$ -l m_mem_free=50G

Rscript --no-save habermann_Macrophage_esvd.R
