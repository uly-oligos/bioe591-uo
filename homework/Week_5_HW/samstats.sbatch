#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-bioe-591-genomics        #specify the account to use
#SBATCH --job-name=samstats                             # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=1              # number of cores to allocate
#SBATCH --time=0-00:30:00                 # Maximum job run time
#SBATCH --output=samstats.out
#SBATCH --error=samstats.err

source ~/.bashrc

cd /home/t62f198/bioe-591-genomics/students/ulysses_oles/Indexing

conda activate samtools

samtools flagstat Diglossa_baritula_77283.sorted.bam > Diglossa_baritula_77283_alignment_stats.txt
samtools depth Diglossa_baritula_77283.sorted.bam | awk '{sum+=$3} END {print sum/NR}' > Diglossa_baritula_77283_depth.txt
