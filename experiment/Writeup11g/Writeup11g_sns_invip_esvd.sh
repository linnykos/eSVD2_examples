#!/bin/bash
#$ -N sns_invip
#$ -j y
#$ -o ../../../../out/Writeup11g/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11g_sns_invip_esvd.R