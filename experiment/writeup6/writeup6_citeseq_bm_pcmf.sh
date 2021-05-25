#!/bin/bash
#$ -N bm_pcmf
#$ -j y
#$ -o ../../../../out/writeup6/qsub/
#$ -pe openmpi 4
#$ -l m_mem_free=12G

Rscript --no-save writeup6_citeseq_bm_pcmf.R
