#!/bin/sh
#PBS -l nodes=1:ppn=1

module purge
module load R/4.0.3

cd $PBS_O_WORKDIR

R --vanilla --file=r2jags-test-model.R