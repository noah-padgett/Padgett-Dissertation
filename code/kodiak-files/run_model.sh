#!/bin/sh
#PBS -l nodes=1:ppn=2

module purge
module load R/3.5.1

cd $PBS_O_WORKDIR

R --vanilla --file=run_model_par.R
