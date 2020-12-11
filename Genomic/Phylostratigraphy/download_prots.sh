#!/bin/bash -l

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=20:0:0

# Request 1 gigabyte of RAM per process.
#$ -l mem=10G

# Request scratch space.
#$ -l tmpfs=15G

# Set the name of the job.
#$ -N download_prots

# Set the working directory
#$ -wd /lustre/home/ucfaata/Scratch/lf_2020/phylostratigraphy

# Load the R module
module unload compilers
module unload mpi
module load r/recommended

# Define path for R libraries
export R_LIBS=/lustre/home/ucfaata/opt/myRlibs

# Run R job
time(Rscript download_prots.R > prots.out)