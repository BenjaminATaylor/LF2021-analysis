  GNU nano 2.4.2                                     File: orthofinder.sh                                                                                  

#!/bin/bash -l

# Request wallclock time (format hours:minutes:seconds).
#$ -l h_rt=5:0:0

# Request RAM.
#$ -l mem=10G

# Request TMPDIR space (default is 10 GB).
#$ -l tmpfs=10G

# Select the number of threads.
#$ -pe mpi 8

# Set the name of the task.
#$ -N Orthofinder_multispecies

# Set working directory
#$ -wd /home/ucfaata/scratch/lf_2020/GO/multispecies

# Set error path
#$ -e /home/ucfaata/scratch/lf_2020/GO/multispecies

# Install gffread onto the cluster
module load python/2.7.9
module load boost/1_54_0/gnu-4.9.2
module load samtools/0.1.19
module load cufflinks/2.2.1

# Install Orthofinder dependencies
module add mcl/14-137
module add muscle/3.8.31

#  Run Orthofinder
# -t: number of threads
# -f: directory of FASTA files for the species we want to analyse
# -A: program to use for MSA (multiple sequence alignment) inference
# -S: program to use for all-versus-all searches instead of BLAST
# -M: method for gene tree inference
time(orthofinder -t 8 -f /home/ucfaata/scratch/lf_2020/GO/multispecies/species_fastas -A muscle -S diamond -M msa)
