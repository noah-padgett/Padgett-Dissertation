#!/bin/sh

module purge
module load R/3.2.0

cd $PBS_O_WORKDIR

R --vanilla --file=howdy.R
