#!/bin/bash
#$ -N sns_astpp
#$ -j y
#$ -o ../../../../out/Writeup13/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_astpp_esvd.R
