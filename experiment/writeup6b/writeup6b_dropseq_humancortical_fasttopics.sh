#!/bin/bash
#$ -N cortical_fasttopics
#$ -j y
#$ -o ../../../../out/writeup6b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup6b_dropseq_humancortical_fasttopics.R
