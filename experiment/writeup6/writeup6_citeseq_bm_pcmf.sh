#!/bin/bash
#$ -N pcmf
#$ -j y
#$ -o ../../../../out/writeup6/qsub/

Rscript --no-save writeup6_citeseq_bm_pcmf.R
