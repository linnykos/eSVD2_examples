#!/bin/bash
#$ -N hab_T_esvd
#$ -j y
#$ -o ../../../../out/Writeup13/qsub/
#$ -l m_mem_free=10G

Rscript --no-save habermann_T_esvd.R
