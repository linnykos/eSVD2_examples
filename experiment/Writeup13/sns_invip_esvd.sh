#!/bin/bash
#$ -N sns_invip
#$ -j y
#$ -o ../../../../out/Writeup13/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_invip_esvd.R
