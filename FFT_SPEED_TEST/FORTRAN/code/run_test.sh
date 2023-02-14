#!/bin/bash
#PBS -A NCGD0011
#PBS -N QG
#PBS -j oe
#PBS -q premium
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=36:ompthreads=36

module purge
module load intel
module load mkl
export OMP_NUM_THREADS=36
export OMP_STACKSIZE=10G
ulimit -s unlimited
ifort -o test.x mkl_dfti.f90 test.f90 -assume byterecl -mkl=parallel -qopenmp -parallel -O3 -xHost -mcmodel medium -shared-intel

### Run OpenMP program
./test.x
