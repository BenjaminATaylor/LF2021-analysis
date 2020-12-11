#!/bin/bash -l

# Request wallclock time (format hours:minutes:seconds).
#$ -l h_rt=5:0:0

# Request RAM.
#$ -l mem=10G

# Request TMPDIR space (default is 10 GB).
#$ -l tmpfs=10G

# Select the number of threads.
#$ -pe mpi 1

# Set the name of the task.
#$ -N BUSCO_LF

# Set working directory
#$ -wd /home/btaylor/LF_analyses/BUSCO

# Set error path
#$ -e /home/btaylor/LF_analyses/BUSCO

source /share/apps/source_files/conda.source
conda activate busco

time(busco -f -m genome -i /home/btaylor/data/Liostenogaster_flavolineata.faa \
 -o BUSCO_LF -l hymenoptera_odb10)
