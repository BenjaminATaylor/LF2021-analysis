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
#$ -N BLAST_dmel_vs_lf

# Set working directory
#$ -wd /home/ucfaata/scratch/lf_2020/orthology/RBH_BLAST

# Set error path
#$ -e /home/ucfaata/scratch/lf_2020/orthology/RBH_BLAST

#load BLAST+
module load blast+/2.2.30/intel-2015-update2

# create BLAST database from Pdom FASTA
makeblastdb -in species_fastas/Liostenogaster_flavolineata.faa -dbtype 'prot' -out databases/LF_db

# run BLAST search
# -evalue: equivalent to significance value cutoff (lower values = more stringent matching)
# -outfmt: out put format (6 = tabular)
blastp -num_threads 8 -evalue 1 -outfmt 6 -query species_fastas/Drosophila_melanogaster.faa -db databases/LF_db > Dmel_to_LF