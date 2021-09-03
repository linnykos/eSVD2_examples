#!/bin/bash
#$ -N pancreas_nb2
#$ -j y
#$ -o ../../../../out/writeup8/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writeup6_10x_mousepancreas_esvd_nb2.R
