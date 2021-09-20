#!/bin/bash
#$ -N pancreas_esvd
#$ -j y
#$ -o ../../../../out/writeup8b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8b_10x_mousepancreas_esvd.R
