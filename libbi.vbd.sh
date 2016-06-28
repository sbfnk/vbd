#!/bin/bash
#$ -N vbd
#$ -V -cwd
#$ -M sebastian.funk@lshtm.ac.uk
#$ -m beas
#$ -l mem_free=2000M,h_vmem=8000M
#$ -pe smp 4
#$ -q parallel.q
#$ -R y

seed=$RANDOM

echo Seed: $seed

## erlang_v=$(echo "(${SGE_TASK_ID} - 1) / 10 + 1" | bc)
## erlang_h=$(echo "(${SGE_TASK_ID} - 1) % 10 + 1" | bc)

##Rscript ~/code/vbd/R/dengue_zika_mcmc.r -n 10000000 -p 10000 -t 4 -e $seed -a -i 1000 -f -r -l 
##Rscript ~/code/vbd/R/dengue_zika_mcmc.r -n 50000 -p 10000 -a -i  1000 -r -l -k -f -v -e 11724 -c 32
Rscript ~/code/vbd/R/dengue_zika_mcmc.r --nsamples 100000 --pre_samples 10000 --seed $seed --sample-prior --sample-observations --keep --force --thin 100 --sero --beta

