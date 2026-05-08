library(tidyverse)

# read sample names and extract
fam <- read_table("mudpuppies.int.fam", 
                  col_names = FALSE,
                  show_col_types = FALSE
)
samples <- fam$X2   # individual IDs are usually column 2

# choose your K value (will be 2 for this demo)
K <- 2

# read Q matrix
q <- read_table("mudpuppies.int.2.Q",
                col_names = FALSE,
                show_col_types = FALSE
)

# name the ancestry columns
colnames(q) <- paste0("Cluster", 1:K)

# combine with sample names
q_df <- q %>%
  mutate(sample = samples) %>%
  relocate(sample)

# convert to long format for ggplot
q_long <- q_df %>%
  pivot_longer(
    cols = starts_with("Cluster"),
    names_to = "cluster",
    values_to = "ancestry"
  ) %>%
  mutate(sample = factor(sample, levels = samples))

# plot
ggplot(q_long, aes(x = sample, y = ancestry, fill = cluster)) +
  geom_col(width = 1, color="white") +
  theme_bw() +
  labs(x = "Individual", y = "Ancestry proportion") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

cv_df <- tibble(file = "admixture.out") %>%
  mutate(
    text = map_chr(file, read_file),
    cv_line = str_extract(text, "CV error \\(K=\\d+\\):\\s*[-0-9.eE]+"),
    K  = str_match(cv_line, "CV error \\(K=(\\d+)\\)")[,2] |> as.integer(),
    CV = str_match(cv_line, ":\\s*([-0-9.eE]+)")[,2] |> as.numeric()
  ) %>%
  select(file, K, CV) %>%
  arrange(K)
cv_df

dna <- vcfR2DNAbin(vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
mudpuppies_genind <- DNAbin2genind(dna)
mudpuppies_genind

mudpuppies_genind_scaled <- scaleGen(mudpuppies_genind,NA.method="mean",scale=F)
mudpuppies_pca <- prcomp(mudpuppies_genind_scaled, center=F,scale=F)

screeplot(mudpuppies_pca)

pc <- data.frame(mudpuppies_pca$x[,1:3])
pc$sample <- rownames(pc)

ggplot(data=pc,aes(x=PC1,y=PC2))+
  geom_text(aes(label=sample))

grp <- find.clusters(mudpuppies_genind, n.pca = 50, n.clust = 2)
grp

pc$cluster <- grp$grp
ggplot(pc, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2) +
  geom_text(aes(label = sample), vjust = -0.5, size = 3) 

dapc2 <- dapc(mudpuppies_genind, pop = grp$grp, n.pca = 50, n.da = 2)
dapc2

q <- as.data.frame(dapc1$posterior)
q$sample <- rownames(q)
q_long <- q |>
  pivot_longer(
    cols = -sample,
    names_to = "cluster",
    values_to = "ancestry"
  )
q_long

ggplot(q_long, aes(x = sample, y = ancestry, fill = cluster)) +
  geom_col(width = 1, color = "white") +
  theme_bw() +
  labs(x = "Individual", y = "Assignment probability") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# identify clusters
grp_auto <- find.clusters(
  mudpuppies_genind,
  n.pca = 50,
  choose.n.clust = FALSE, 
  max.n.clust = 10,
  stat = "BIC"
)

# print BIC values
grp_auto$Kstat

# plot!
plot(
  1:length(grp_auto$Kstat),
  grp_auto$Kstat,
  type = "b",
  xlab = "K",
  ylab = "BIC"
)


####Homework

K <- 3

q <- read_table("mudpuppies.int.2.Q",
                col_names = FALSE,
                show_col_types = FALSE
)

# name the ancestry columns
colnames(q) <- paste0("Cluster", 1:K)

# combine with sample names
q_df <- q %>%
  mutate(sample = samples) %>%
  relocate(sample)

# convert to long format for ggplot
q_long <- q_df %>%
  pivot_longer(
    cols = starts_with("Cluster"),
    names_to = "cluster",
    values_to = "ancestry"
  ) %>%
  mutate(sample = factor(sample, levels = samples))

# plot
ggplot(q_long, aes(x = sample, y = ancestry, fill = cluster)) +
  geom_col(width = 1, color="white") +
  theme_bw() +
  labs(x = "Individual", y = "Ancestry proportion") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

cv_df <- tibble(file = "admixture.out") %>%
  mutate(
    text = map_chr(file, read_file),
    cv_line = str_extract(text, "CV error \\(K=\\d+\\):\\s*[-0-9.eE]+"),
    K  = str_match(cv_line, "CV error \\(K=(\\d+)\\)")[,2] |> as.integer(),
    CV = str_match(cv_line, ":\\s*([-0-9.eE]+)")[,2] |> as.numeric()
  ) %>%
  select(file, K, CV) %>%
  arrange(K)
cv_df

library(adegenet)
library(vcfR)

vcf <- read.vcfR(file = "mudpuppies.vcf", verbose = TRUE)

dna <- vcfR2DNAbin(vcf, unphased_as_NA = F, consensus = T, extract.haps = F)
mudpuppies_genind <- DNAbin2genind(dna)
mudpuppies_genind

mudpuppies_genind_scaled <- scaleGen(mudpuppies_genind,NA.method="mean",scale=F)
mudpuppies_pca <- prcomp(mudpuppies_genind_scaled, center=F,scale=F)

screeplot(mudpuppies_pca)

pc <- data.frame(mudpuppies_pca$x[,1:3])
pc$sample <- rownames(pc)

ggplot(data=pc,aes(x=PC1,y=PC2))+
  geom_text(aes(label=sample))

grp <- find.clusters(mudpuppies_genind, n.pca = 50, n.clust = 2)
grp

pc$cluster <- grp$grp
ggplot(pc, aes(x = PC1, y = PC2, color = cluster)) +
  geom_point(size = 2) +
  geom_text(aes(label = sample), vjust = -0.5, size = 3) 

q <- as.data.frame(dapc1$posterior)
q$sample <- rownames(q)
q_long <- q |>
  pivot_longer(
    cols = -sample,
    names_to = "cluster",
    values_to = "ancestry"
  )
q_long

ggplot(q_long, aes(x = sample, y = ancestry, fill = cluster)) +
  geom_col(width = 1, color = "white") +
  theme_bw() +
  labs(x = "Individual", y = "Assignment probability") +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

# identify clusters
grp_auto <- find.clusters(
  mudpuppies_genind,
  n.pca = 50,
  choose.n.clust = FALSE, 
  max.n.clust = 10,
  stat = "BIC"
)

# print BIC values
grp_auto$Kstat

# plot!
plot(
  1:length(grp_auto$Kstat),
  grp_auto$Kstat,
  type = "b",
  xlab = "K",
  ylab = "BIC"
)