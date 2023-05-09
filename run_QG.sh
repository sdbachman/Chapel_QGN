#!/bin/bash
#PBS -A NCGD0011
#PBS -N ChapQG
#PBS -j oe
#PBS -q premium
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=36:mpiprocs=36:mem=109GB
#PBS -o /glade/scratch/bachman/Ian_QG/QGN_Beta/TEST/CHPL

module purge
source ~/.chapel_QG_64

#./QGN_Driver --restart=false --Nt=10000 --Q_DIAG_FREQ=50 --QG_Leith_coeff=500.0
./QGN_Driver --restart=false --Nt=10000 --nz=3 --background_file=background_eady.nc1 --Q_DIAG_FREQ=50 --QG_Leith_coeff=500.0
