#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-bioe-591-genomics        #specify the account to use
#SBATCH --job-name=HW7het.sbatch                             # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=1              # number of cores to allocate
#SBATCH --time=0-02:00:00                 # Maximum job run time
#SBATCH --output=HW7het.out
#SBATCH --error=HW7het.err

source ~/.bashrc

cd /home/t62f198/bioe-591-genomics/students/ulysses_oles/Indexing

conda activate vcftools

vcftools --vcf Diglossa_baritula_multisample_cleaned.recode.vcf --het
