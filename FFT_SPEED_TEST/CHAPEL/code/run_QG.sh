#!/bin/bash
#PBS -A NCGD0011
#PBS -N QG
#PBS -j oe
##PBS -k eod
##PBS -m abe
##PBS -M email_address
#PBS -q regular
#PBS -l walltime=12:00:00
#PBS -l select=1:ncpus=36

module purge
source ~/.chapel_pbs-gasnet

REPL1
