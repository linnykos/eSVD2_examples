#!/bin/bash
#$ -N sns_invip2
#$ -j y
#$ -o ../../../../out/Writeup11e/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11e_sns_invip_esvd2.R
