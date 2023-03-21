#!/bin/bash
#$ -N sns_insst
#$ -j y
#$ -o ../../../../out/Writeup13/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_insst_esvd.R
