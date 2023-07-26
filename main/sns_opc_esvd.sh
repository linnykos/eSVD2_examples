#!/bin/bash
#$ -N sns_opc_esvd
#$ -j y
#$ -o ../../../out/main/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_opc_esvd.R
