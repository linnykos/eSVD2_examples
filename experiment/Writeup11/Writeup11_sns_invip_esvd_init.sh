#!/bin/bash
#$ -N sns_invip_esvd_init
#$ -j y
#$ -o ../../../../out/Writeup11/qsub/
#$ -l m_mem_free=50G

Rscript --no-save Writeup11_sns_invip_esvd_init.R