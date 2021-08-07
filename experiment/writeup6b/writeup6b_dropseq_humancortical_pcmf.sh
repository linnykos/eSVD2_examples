#!/bin/bash
#$ -N cortical_pcmf
#$ -j y
#$ -o ../../../../out/writeup6b/qsub/
#$ -pe openmpi 4
#$ -l m_mem_free=50G

Rscript --no-save writeup6b_dropseq_humancortical_pcmf.R
