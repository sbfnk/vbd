#!/bin/bash
#$ -N vbd
#$ -V -cwd
#$ -M sebastian.funk@lshtm.ac.uk
#$ -l mem_free=2000M,h_vmem=8000M
#$ -pe smp 8
#$ -q parallel.q
#$ -R y
#$ -t 1-100

seed=$RANDOM

echo Seed: $seed

Rscript ~/code/vbd/R/dengue_zika_mcmc.r --nsamples 1000000 --pre_samples 1000 --sample-prior --sample-observations  --force --thin 1000 --patch --poisson --seed $seed --parallel-number ${SGE_TASK_ID} --sero

