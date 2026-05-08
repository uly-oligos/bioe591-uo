library(ape)
library(readr)

# read tree
tree <- read.tree("lemurs.snps.min4.phy.contree")
tree

tree$tip.label

# read metadata
meta <- read_tsv("lemur_metadata.txt")

# map species to sample
map <- setNames(meta$species, meta$ID_long)

# overwrite tip labels
tree$tip.label <- unname(map[tree$tip.label])

plot(tree)

# root using an outgroup (replace with your sample name)
tree_rooted <- root(tree, outgroup = "murinus", resolve.root = TRUE)

# plot
plot(tree_rooted)

# add node labels
bs <- as.numeric(tree$node.label)
nodelabels(ifelse(bs >= 70, bs, ""), cex = 0.7, frame = "n")

png("lemur_tree.png", width = 1000, height = 800)
plot(tree_rooted)
nodelabels(tree$node.label, cex = 0.7, frame = "n")
dev.off()

#My analysis does not support the lumping of M. mittermeieri and M. lehilahytsara. My tree does, however
#support the idea that species 3 is distinct