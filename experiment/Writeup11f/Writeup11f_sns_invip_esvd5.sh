#!/bin/bash
#$ -N sns_invip5
#$ -j y
#$ -o ../../../../out/Writeup11f/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11f_sns_invip_esvd5.R
