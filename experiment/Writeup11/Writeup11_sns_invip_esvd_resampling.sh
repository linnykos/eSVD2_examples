#!/bin/bash
#$ -N sns_invip_esvd_resample
#$ -j y
#$ -o ../../../../out/Writeup11/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11_sns_invip_esvd_resampling.R