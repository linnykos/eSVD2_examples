#!/bin/bash
#$ -N habermann_preprocess_step1
#$ -j y
#$ -o ../../../../out/Writeup11/qsub/
#$ -l m_mem_free=100G

Rscript --no-save Writeup11_habermann_preprocess_step1.R
