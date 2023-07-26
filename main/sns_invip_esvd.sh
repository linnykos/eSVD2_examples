#!/bin/bash
#$ -N sns_invip_esvd
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_invip_esvd.R
