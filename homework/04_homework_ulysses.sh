#!/bin/bash
##
## example-array.slurm.sh: submit an array of jobs with a varying parameter
##
## Lines starting with #SBATCH are read by Slurm. Lines starting with ## are comments.
## All other lines are read by the shell.
##
#SBATCH --account=priority-bioe-591-genomics        #specify the account to use
#SBATCH --job-name=fastp_loop                             # job name
#SBATCH --partition=priority              # queue partition to run the job in
#SBATCH --nodes=1                       # number of nodes to allocate
#SBATCH --ntasks-per-node=1             # number of descrete tasks - keep at one except for MPI
#SBATCH --cpus-per-task=1              # number of cores to allocate
#SBATCH --time=0-00:30:00                 # Maximum job run time
#SBATCH --output=fastp_loop.out
#SBATCH --error=fastp_loop.err

source ~/.bashrc

cd /home/t62f198/bioe-591-genomics/course-materials/data/raw_reads/Diglossa_baritula/

# activate your fastp environment
conda activate fastp

for i in `ls -1 *R1_001.fastq.gz | sed 's/R1_001.fastq.gz//'`
do
fastp -i $i\R1_001.fastq.gz \
-o ~/bioe-591-genomics/students/ulysses_oles/trimmed_reads/$i\R1_001.trimmed.fastq.gz \
-I $i\R2_001.fastq.gz \
-O ~/bioe-591-genomics/students/ulysses_oles/trimmed_reads/$i\R2_001.fastq.gz \
-h ~/bioe-591-genomics/students/ulysses_oles/trimmed_reads/$i\.fastp.html \
-j ~/bioe-591-genomics/students/ulysses_oles/trimmed_reads/$i\.fastp.json
done
