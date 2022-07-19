#!/bin/bash
#$ -N regev_preprocess
#$ -j y
#$ -o ../../../../out/Writeup11h/qsub/
#$ -l m_mem_free=100G

Rscript --no-save regev_preprocess_celltypes.R
