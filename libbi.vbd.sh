#!/bin/bash
#$ -N vbd
#$ -V -cwd
#$ -M sebastian.funk@lshtm.ac.uk
#$ -l mem_free=2000M,h_vmem=8000M
#$ -pe smp 4
#$ -q parallel.q
#$ -R y
#$ -t 1-100

seed=$RANDOM

echo Seed: $seed

##Rscript ~/code/vbd/R/dengue_zika_mcmc.r -n 10000000 -p 10000 -t 4 -e $seed -a -i 1000 -f -r -l 
##Rscript ~/code/vbd/R/dengue_zika_mcmc.r -n 50000 -p 10000 -a -i  1000 -r -l -k -f -v -e 11724 -c 32
Rscript ~/code/vbd/R/dengue_zika_mcmc.r --nsamples 100000 --pre_samples 10000 --seed $seed --sample-prior --sample-observations --keep --force --thin 100 --sero --beta
Rscript ~/code/vbd/R/dengue_zika_mcmc.r --nsamples 1000000 --pre_samples 1000 --sample-prior --sample-observations  --force --thin 1000 --patch --poisson --seed $seed --parallel-number ${SGE_TASK_ID}

