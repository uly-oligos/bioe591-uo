#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-bioe-591-genomics        #specify the account to use
#SBATCH --job-name=HW7.sbatch                             # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=1              # number of cores to allocate
#SBATCH --time=0-02:00:00                 # Maximum job run time
#SBATCH --output=HW7.out
#SBATCH --error=HW7.err

source ~/.bashrc

cd /home/t62f198/bioe-591-genomics/students/ulysses_oles/Indexing

conda activate vcftools

vcftools --gzvcf Diglossa_baritula_multi_sample.vcf.gz --recode --out Diglossa_baritula_multisample_cleaned --remove-indels --minQ 40 --mac 1 --max-missing-count .1 --min-meanDP 1 --thin 50
