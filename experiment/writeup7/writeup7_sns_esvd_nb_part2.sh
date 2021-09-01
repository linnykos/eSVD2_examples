#!/bin/bash
#$ -N sns_esvd_nb_part2
#$ -j y
#$ -o ../../../../out/writeup7/qsub/
#$ -l m_mem_free=100G

Rscript --no-save writeup7_sns_esvd_nb_part2.R
