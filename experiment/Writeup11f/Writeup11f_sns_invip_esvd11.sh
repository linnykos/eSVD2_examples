#!/bin/bash
#$ -N sns_invip11
#$ -j y
#$ -o ../../../../out/Writeup11f/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11f_sns_invip_esvd11.R
