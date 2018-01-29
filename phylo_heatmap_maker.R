# Load packages
library("ape")
library(gplots)

# All input for the whole script
# Tree in newich format
tree_file <- "C:/Users/Arne/Desktop/Thesis/phylogeny/complete_Streptomyces_tree.nwk"
# Tree_outgroup needs to be an exact node name needed to root the tree
tree_outgroup <- "Salinispora_arenicola_CNS-205_CP000850.1_CP000850"
# Matrix for the heatmap
absence_presence_file <- "C:/Users/Arne/Desktop/Thesis/phylogeny/heatmap_input_c80_other.csv"
# Minimum number of species a BGC must be in to be in the more abundant plot later
# after the plot with all BGCs
min_amount_species <- 5

#Prepare the phylogenetic tree
mytree <- read.tree(tree_file)
mytree_rooted <- root(mytree, outgroup = tree_outgroup,resolve.root = T)
mytree_rooted$edge.length[which(mytree_rooted$edge.length == 0)] <- 0.00001
mytree_um <- chronopl(mytree_rooted,lambda = 0.1, tol = 0)
mytree_dendro <- as.dendrogram(as.hclust.phylo(mytree_um))

#force row order so that it matches the order of leafs in rep_tree_d
clade_order <- order.dendrogram(mytree_dendro)
clade_name <- labels(mytree_dendro)
clade_position <- data.frame(clade_name, clade_order)
clade_position <- clade_position[order(clade_position$clade_order),]

#Load in the absence presence matrix
my_absence_presence <- read.csv(absence_presence_file, header=TRUE)
names <- as.character(my_absence_presence$species[1:100])
my_absence_presence <- my_absence_presence[1:100,2:ncol(my_absence_presence)]
# colnames(my_absence_presence) <- strtrim(colnames(my_absence_presence), 10)
my_absence_presence <- as.matrix(my_absence_presence)
rownames(my_absence_presence) <- names
new_order <- match(clade_position$clade_name, row.names(my_absence_presence))
ordered_absence_presence <- my_absence_presence[new_order,]
# 
# #plot the heatmap
# heatmap.2(ordered_absence_presence, Rowv=mytree_dendro, Colv=F, dendrogram='row',col =
#             colorRampPalette(c("white","red"))(2),
#           sepwidth=c(0.01,0.02),sepcolor="black",colsep=1:ncol(my_absence_presence),rowsep=1:nrow(my_absence_presence),
#           key=FALSE,trace="none",
#           cexRow=.30,cexCol=.5,srtCol=90,
#           margins=c(10,10),
#           main="BGC presence and absence")

# make a subset of columns that appear in more than "x" species
interesting_columns <- apply(my_absence_presence, 2, sum)
interesting_columns_matrix <- ordered_absence_presence[,interesting_columns >= min_amount_species]

#plot the heatmap with more abundannt columns
heatmap.2(interesting_columns_matrix, Rowv=mytree_dendro, Colv=F, dendrogram='row',col = colorRampPalette(c("white","red"))(2),
          sepwidth=c(0.01,0.01),sepcolor="black",rowsep=1:nrow(my_absence_presence),
          key=FALSE,trace="none",
          cexRow=.5,cexCol=1,srtCol=30,
          margins=c(10,10), lwid = c(.5,4), lhei = c(0.5,6),
          main="More abundant BGCs")
