#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-bioe-591-genomics        #specify the account to use
#SBATCH --job-name=alignagain                             # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=1              # number of cores to allocate
#SBATCH --time=0-00:30:00                 # Maximum job run time
#SBATCH --output=alignagain.out
#SBATCH --error=alignagain.err

source ~/.bashrc

cd /home/t62f198/bioe-591-genomics/students/ulysses_oles/Indexing

conda activate align

bwa mem -t 4 hemoglobin_references.fasta /home/t62f198/bioe-591-genomics/students/ulysses_oles/trimmed_reads/Diglossa_baritula_77283_R1_001.trimmed.fastq.gz /home/t62f198/bioe-591-genomics/students/ulysses_oles/trimmed_reads/Diglossa_baritula_77283_R2_001.trimmed.fastq.gz > Diglossa_baritula_77283.sam
