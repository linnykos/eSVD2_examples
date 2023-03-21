#!/bin/bash
#$ -N sns_opc
#$ -j y
#$ -o ../../../../out/Writeup13/qsub/
#$ -l m_mem_free=15G

Rscript --no-save sns_opc_esvd.R
