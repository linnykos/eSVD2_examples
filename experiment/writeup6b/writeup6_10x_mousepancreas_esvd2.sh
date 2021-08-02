#!/bin/bash
#$ -N pancreas_esvd2
#$ -j y
#$ -o ../../../../out/writeup6b/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup6_10x_mousepancreas_esvd2.R
