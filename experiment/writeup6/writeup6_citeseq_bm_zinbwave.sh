#!/bin/bash
#$ -N bm_zinbwave
#$ -j y
#$ -o ../../../../out/writeup6/qsub/
#$ -pe openmpi 4

Rscript --no-save writeup6_citeseq_bm_zinbwave.R
