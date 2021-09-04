#!/bin/bash
#$ -N cortical_esvd
#$ -j y
#$ -o ../../../../out/writeup8/qsub/
#$ -l m_mem_free=50G

Rscript --no-save writupe8_dropseq_humancortical_esvd_nb.R
