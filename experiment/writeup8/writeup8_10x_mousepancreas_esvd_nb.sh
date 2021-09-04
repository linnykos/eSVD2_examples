#!/bin/bash
#$ -N pancreas_nb
#$ -j y
#$ -o ../../../../out/writeup8/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup8_10x_mousepancreas_esvd_nb.R
