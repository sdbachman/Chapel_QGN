#!/bin/bash
#PBS -A NCGD0011
#PBS -N ChapQG
#PBS -j oe
#PBS -q premium
#PBS -l walltime=02:00:00
#PBS -l select=1:ncpus=36:mpiprocs=36:mem=109GB
#PBS -o /glade/scratch/bachman/Ian_QG/QGN_Beta/TEST/CHPL

module purge
source ~/.chapel_QG_64

./QGN_Driver --restart=true --Nt=100 --Q_DIAG_FREQ=50
