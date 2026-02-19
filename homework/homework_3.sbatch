#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-bioe-591-genomics        #specify the account to use
#SBATCH --job-name=homework_3                             # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=1              # number of cores to allocate
#SBATCH --time=0-00:10:00                 # Maximum job run time
#SBATCH --output=seqtk_hw3-%j.out
#SBATCH --error=seqtk_hw3-%j.err

## Run 'man sbatch' for more information on the options above.

source ~/.bashrc

conda activate demo_new

seqtk sample -s100 ~/bioe-591-genomics/course-materials/data/intro/catamenia_analis.fastq.gz 1000 > /home/t62f198/bioe-591-genomics/students/ulysses_oles/subsampled_hw3.fq
