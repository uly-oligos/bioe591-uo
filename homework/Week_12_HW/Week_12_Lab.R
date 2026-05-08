install.packages("pcadapt")
install.packages("patchwork")
install.packages("viridis")
library(pcadapt)
library(vcfR)
library(adegenet)
library(ggplot2)
library(dplyr)
library(readr)
library(patchwork)
library(viridis)

meta <- read_csv("metadata.csv")
x <- read.pcadapt("salvelinus.vcf", type = "vcf")

pca <- pcadapt(x, K = 20)
plot(pca, option = "screeplot")

pca_4 <- pcadapt(x, K = 4)

alpha <- 0.05
padj <- p.adjust(pca_4$pvalues, method = "bonferroni")
outlier_idx <- which(padj < alpha)
neutral_idx <- which(padj >= alpha)
cat("Outlier loci:", length(outlier_idx), "\n")

cat("Neutral loci:", length(neutral_idx), "\n")

plot(pca_4)

vcf <- read.vcfR("salvelinus.vcf", verbose = FALSE)
vcf_neutral <- vcf[neutral_idx, ]
vcf_outlier <- vcf[outlier_idx, ]

dna_neutral <- vcfR2DNAbin(vcf_neutral, unphased_as_NA = FALSE, consensus = TRUE, extract.haps = FALSE)
gi_neutral <- DNAbin2genind(dna_neutral)
dna_outlier <- vcfR2DNAbin(vcf_outlier, unphased_as_NA = FALSE, consensus = TRUE, extract.haps = FALSE)
gi_outlier <- DNAbin2genind(dna_outlier)

gi_neutral_scaled <- scaleGen(gi_neutral, NA.method = "mean", scale = FALSE)
pca_neutral <- prcomp(gi_neutral_scaled, center = FALSE, scale. = FALSE)
gi_outlier_scaled <- scaleGen(gi_outlier, NA.method = "mean", scale = FALSE)
pca_outlier <- prcomp(gi_outlier_scaled, center = FALSE, scale. = FALSE)

# turn lists into dataframes
neutral_df <- data.frame(pca_neutral$x[, 1:2]) %>%
  mutate(sample = rownames(.))
outlier_df <- data.frame(pca_outlier$x[, 1:2]) %>%
  mutate(sample = rownames(.))

# join metadata
neutral_df <- left_join(neutral_df, meta, by = c("sample" = "Sample.ID"))
outlier_df <- left_join(outlier_df, meta, by = c("sample" = "Sample.ID"))

# shared axis limits
xlims <- range(c(outlier_df$PC1, neutral_df$PC1), na.rm = TRUE)
ylims <- range(c(outlier_df$PC2, neutral_df$PC2), na.rm = TRUE)

# percent variance explained
neutral_var <- 100 * summary(pca_neutral)$importance[2, ]
outlier_var <- 100 * summary(pca_outlier)$importance[2, ]

#PCA plot comparison for both sets of loci relative to latitude
p1 <- ggplot(neutral_df, aes(x = PC1, y = PC2, color=Site.Latitude)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_viridis() +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(
    title = "PCA: neutral loci",
    x = paste0("PC1 (", round(neutral_var[1], 1), "%)"),
    y = paste0("PC2 (", round(neutral_var[2], 1), "%)")
  ) 

p2 <- ggplot(outlier_df, aes(x = PC1, y = PC2, color=Site.Latitude)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_viridis() +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(
    title = "PCA: outlier loci only",
    x = paste0("PC1 (", round(outlier_var[1], 1), "%)"),
    y = paste0("PC2 (", round(outlier_var[2], 1), "%)")
  )

p1 + p2

### The neutral loci in the dataset seem most strongly to differentiate along PC2, but this particular
# axis does not explain a very large percentage of the differentiation. Looking at the outlier loci
# specifically, there appears to be a much stronger correlation between PC1 and the site latitude.This
# suggests to me that differentiation along latitude is present, but not the strongest driver in the neutral
# loci, while there is a very strong connection between latitude and variation in the outlier loci.


#This plots against longitude
p1 <- ggplot(neutral_df, aes(x = PC1, y = PC2, color=Site.Longitude)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_viridis() +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(
    title = "PCA: neutral loci",
    x = paste0("PC1 (", round(neutral_var[1], 1), "%)"),
    y = paste0("PC2 (", round(neutral_var[2], 1), "%)")
  ) 

p2 <- ggplot(outlier_df, aes(x = PC1, y = PC2, color=Site.Longitude)) +
  geom_point(size = 2, alpha = 0.8) +
  scale_color_viridis() +
  coord_cartesian(xlim = xlims, ylim = ylims) +
  theme_classic() +
  theme(legend.position = "bottom") +
  labs(
    title = "PCA: outlier loci only",
    x = paste0("PC1 (", round(outlier_var[1], 1), "%)"),
    y = paste0("PC2 (", round(outlier_var[2], 1), "%)")
  )

p1 + p2

##### Comparing these two PCAs showing longitude data, I feel that they both suggest that longitude
# is not a very strong driver of differentiation in either set of loci. There is potentially a weak
# association between longitude and differences in neutral loci, and no obvious pattern of 
# correlation between the two in outlier loci.
