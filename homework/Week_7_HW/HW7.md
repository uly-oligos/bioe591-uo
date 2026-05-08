--remove-indels is used to remove insertions and deletions from the data that can complicate downstream analyses
-minQ 40 filters out reads with a lower phred quality score than 40 to help ensure that artifacts are removed from the dataset
--mac .1 removes all loci with a minor allele count of less than .1 
--max-missing-count 1 means I removed all sites with missing data
--min-meanDP 1 only keeps sites that have a minimum average depth of 1
--thin 50 thins the dataset to one variant site randomly every 50 bases

I ran --het for the dataset and ended up with file containing to heterozygosity information. I believe this stems from me overfiltering the dataset massively
