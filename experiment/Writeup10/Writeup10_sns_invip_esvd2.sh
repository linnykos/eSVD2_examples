#!/bin/bash
#$ -N sns_invip_esvd2
#$ -j y
#$ -o ../../../../out/Writeup10/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup10_sns_invip_esvd2.R