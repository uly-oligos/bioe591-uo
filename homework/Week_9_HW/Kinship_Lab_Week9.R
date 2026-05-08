install.packages("related", repos="http://R-Forge.R-project.org")
install.packages("adegenet")
install.packages("vcfR")
install.packages("pegas")

library(related)
library(adegenet)
library(vcfR)
library(pegas)
library(tidyr)
library(reshape2)
library(ggplot2)
library(tidyverse)

vcf <- read.vcfR(file = "Inca_MaxMissing10.recode.vcf", verbose = TRUE)

vcf
vcf@meta
vcf@fix
vcf@gt

genind_obj <- vcfR2genind(vcf)

genind_obj
head(genind_obj@tab)
summary(genind_obj@loc.n.all)

af_summary <- adegenet::summary(genind_obj) # here it is important to specify which package's summary function is used!
af_summary # view object
h_o <- af_summary$Hobs
h_e <- af_summary$Hexp

head(h_o)
head(h_e)
het_df <- data.frame(locus = names(h_o), h_o = h_o, h_e = h_e)
head(het_df)

Fis_per_locus <- 1 - (h_o / h_e)
Fis_per_locus
mean(Fis_per_locus, na.rm = TRUE)

loci_obj <- genind2loci(genind_obj)
hwe_results <- pegas::hw.test(loci_obj, B = 100)
hwe_results

hwe_results %>% as_tibble() %>% filter(Pr.exact<0.05)

gt_filtered <- vcfR::extract.gt(vcf, element = "GT")

# sample ids
sample_ids <- colnames(gt_filtered)

gt_to_alleles <- function(gt_vector) {
  # split "0/1" or "0|1" into two integer alleles, returning a 2-column matrix (samples x 2 alleles)
  allele1 <- integer(length(gt_vector))
  allele2 <- integer(length(gt_vector))
  
  for (i in seq_along(gt_vector)) {
    g <- gt_vector[i]
    if (is.na(g) || g %in% c("./.", ".", "./", "/.")) {
      allele1[i] <- 0
      allele2[i] <- 0
    } else {
      parts <- as.integer(strsplit(g, "[/|]")[[1]])
      allele1[i] <- parts[1] + 1L    # shift: 0->1 (ref), 1->2 (alt)
      allele2[i] <- parts[2] + 1L
    }
  }
  cbind(allele1, allele2)
}

allele_list <- vector("list", nrow(gt_filtered))

for (v in seq_len(nrow(gt_filtered))) {
  allele_list[[v]] <- gt_to_alleles(gt_filtered[v, ])
}

# combine: each element is (n_samples x 2); bind column-wise
allele_matrix <- do.call(cbind, allele_list)

# add individual IDs as the first column
coancestry_input <- data.frame(IndID = sample_ids, allele_matrix,
                               stringsAsFactors = FALSE)

# column names: IndID, L1_a, L1_b, L2_a, L2_b, ...
locus_names <- paste0(rep(paste0("L", seq_len(nrow(gt_filtered))),
                          each = 2),
                      rep(c("_a", "_b"), nrow(gt_filtered)))
colnames(coancestry_input) <- c("IndID", locus_names)

coancestry_input[1:5, 1:7]

kin_results <- related::coancestry(
  genotype.data = coancestry_input,
  wang          = 1,      # 1 = compute; 0 = skip
  dyadml        = 1,
  quellergt     = 1
)

head(kin_results$relatedness)

###Homework

ngsrelate.data <- read.table("Inca_MaxMissing10.ngsrelate.out", header = TRUE)

average_kinship_ngsrelate <- mean(ngsrelate.data$theta)
average_kinship_ngsrelate

average_kinship_related <- mean(kin_results$relatedness$dyadml)
average_kinship_related

###These functions calculate the average kinship coefficients between the two datasets. 
###This comparison shows a higher average coefficient of kinship in the related dataset