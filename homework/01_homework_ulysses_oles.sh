#!/usr/bin/env bash
cd ~/Desktop/test/week_1
mkdir fastq/
mkdir fasta/
mkdir metadata/

mv *.fastq.gz fastq/
mv *fasta fasta/
mv *csv metadata/

ls -1 ~/Desktop/test/week_1/fastq/ | wc -l
ls -1 ~/Desktop/test/week_1/fasta/ | wc -l
ls -1 ~/Desktop/test/week_1/metadata/ | wc -l

echo "DONE!"