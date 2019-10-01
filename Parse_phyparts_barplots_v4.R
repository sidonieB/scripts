###################################################################################################
# This script parses the results of phyparts (from Steven Smith) to obtain a classification of genes
# depending on how they conflict and depending on a user-defined categorization, e.g. plastid/nuclear
# It requires to have run 2 phyparts analyses, with different -s settings (here no setting + -s 70)
# Was used to produce results for Bellot et al. 2019. Scientific Reports
# Probably not ready-to-use by someone else but feel free to try or ask me
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
####################################################################################################


library(RColorBrewer)
library(randomcoloR)
library(utils)
library(ggplot2)#
library(scales)
library(dplyr)
library(ggforce)
library(ape) #
library(phangorn)#
library(ggtree)#
library(ggimage)
library(cowplot)

# Input

setwd("C:/Users/xxxxx/Copies_TreesForPlots_2/")        # set the path to the directory where you want the output
input <- "C:/Users/xxxx/Phyparts/RAxML_A/RAxML_bestTree_A_conc_part_rooted2_gene_summary.txt"        # set the path to the input table, which is the table obtained from the python script
inputKey <- "C:/Users/xxxx/Phyparts/RAxML_A/phyparts_A_BP0.node.key"     # set the path to the input key, which is the *node.key from phyparts
inputTree <- "C:/Users/xxxx/Phyparts/RAxML_A/RAxML_bestTree.A_conc_part_rooted3.tre"     # set the path to the input tree, which is the species tree on which you ran phyparts
#phypartsTreesLS <- "C:/Users/xxxxx/Phyparts/RAxML_A/phyparts_A_BP0.concon.tre"    # set the path to the list of 3 trees outputed by phyparts when using all bipartitions
#phypartsTreesHS <- "C:/Users/xxxxxx/Phyparts/RAxML_A/phyparts_A_BP70.concon.tre"    # set the path to the list of 3 trees outputed by phyparts when collapsing clades with less than a certain BP support
#phyparts_altsLS <- "C:/Users/xxxxxx/Phyparts/RAxML_A/phyparts_A_BP0.hist.alts"     # set the path to the file of alternative bipartitions outputed by phyparts when using all bipartitions
#phyparts_altsHS <- "C:/Users/xxxxxxx/Phyparts/RAxML_A/phyparts_A_BP70.hist.alts"     # set the path to the file of alternative bipartitions outputed by phyparts when collapsing clades with less than a certain BP support
GenesNum <- 277 # number of genes used in the phyparts analyis
Title_pies <- "RAxML, all loci" # title for the tree with pies plot
output_file <- "Figure_RAxML_A_pies_and_bars.pdf" # name of the output file with the plots


# import tree and node key from phypart, 
# and make a table with correspondancy between nodes as in phyparts (NodesP) and nodes as R reads the tree (NodesT)

t <- read.tree(inputTree)
k <- read.table(inputKey, header=F, row.names=NULL, sep=" ", na.strings = c("", "NA"), stringsAsFactors=FALSE)
colnames(k) <- c("NodeP", "Tips")
NodeP <- rep(0, length(k$NodeP))
k2 <- as.data.frame(NodeP)
k2$NodeT <- rep(0, length(k2$NodeP))

for (x in 1:length(k[,1])) {
  k2$NodeP[x] <-  k[x,1]
  tips <- gsub("\\)|\\(", "", k[x,2])
  tips_vec <- strsplit(tips, ",")
  k2$NodeT[x] <- mrca.phylo(t, tips_vec[[1]])
}


# plot tree with node numbered like in phyparts or following R numbering

t$node_numbT <- c("", as.character(k2$NodeT))
t$node_numbP <- c("", as.character(k2$NodeP))

pP <- ggtree(t, ladderize=T) + 
  geom_tiplab(size=2.5) + 
  xlim_tree(0.15) + 
  ggtitle("Species tree, RAxML, phyparts node labels") + 
  geom_nodelab(label=t$node_numbP, size=3, vjust=0, hjust= -0.1)

pT <- ggtree(t, ladderize=T) + 
  geom_tiplab(size=2.5) + 
  xlim_tree(0.15) + 
  ggtitle("Species tree, RAxML, R node labels") + 
  geom_nodelab(label=t$node_numbT, size=3, vjust=0, hjust= -0.1)






##########################################
# Produce barplots for each node with the frequencies of genes that agree or disagree or are neutral regarding this node
# do it with or without taking into account bootstrap support
# indicate if they are involved in a top alternative
# and what is their category
##########################################

# import and format the gene summary table
df <-read.table(input, header=T, row.names=NULL, sep="\t", na.strings = c("", "NA"), stringsAsFactors=FALSE)
as.data.frame(df)
head(df)

# build a "dictionary" (in fact list) of node numbers, with the phyparts numbers as keys and the tree numbers as values
# could have made a hash - anyway...

nodes_dictPkey <- vector(mode="list", length=length(k2$NodeP))
names(nodes_dictPkey) <- as.character(k2$NodeP)
for (z in 1:length(nodes_dictPkey)) {
  nodes_dictPkey[[z]] <- k2$NodeT[z]
}


# add a column to the df with the node number as in the tree

df$nodeT <- rep(0, length(df$node))
for (y in 1:length(df$node)) {
  n <- as.character(df$node[y])    # very important to have it as character or will use the node number as index to find the value in the dictionary
  df$nodeT[y] <- nodes_dictPkey[[n]]
}

# split dataset per node based on nodeT: for each node, we have a subdataset
dfnode <- split(df, df$nodeT) 
length(dfnode)
names(dfnode)


# build a small table with number of rows = category*status*number of nodes
# and column names being nodeT, status, category, Freq

table_rows <- length(unique(df$category)) * length(unique(df$status)) * length(names(dfnode))

df2 <- data.frame(matrix(ncol = 4, nrow = table_rows )) # +1 because of the frequency column of the contingency tables
colnames(df2) <- c("nodeT", "status", "category", "Freq")

l <- c("ABSENT_BIP","AGREE_HS","AGREE_LS","CONFLICT_HS","CONFLICT_LS")

rowTot <- 0

for (s in (1:length(dfnode))) {         # for each subset
  Nuclear <- dfnode[[s]][dfnode[[s]]$category == 'Nuclear',]
  Plastid <- dfnode[[s]][dfnode[[s]]$category == 'Plastid',]
  tabN <- table(factor(Nuclear$status, levels = l)) 
  tabP <- table(factor(Plastid$status, levels = l))
  tabNP <- c(tabN, tabP)
  freq <- rep(0, length(tabNP))
  for (x in 1:length(tabNP)) {freq[x] <- tabNP[[x]]}
  
  df2[(rowTot+1):(rowTot+10),1] <- rep(dfnode[[s]]$nodeT[1], 10)
  df2[(rowTot+1):(rowTot+10),2] <- c(names(tabN), names(tabP))
  df2[(rowTot+1):(rowTot+10),3] <- c(rep("Nuclear", 5), rep("Plastid", 5))
  df2[(rowTot+1):(rowTot+10),4] <- freq
  
  rowTot <- rowTot+10
  
}


p4 <- ggplot(df2, aes(x = status, y = Freq, fill = category)) +
  geom_bar(position = "stack", stat = "identity") +
  facet_wrap( ~ nodeT) + 
  theme(#panel.spacing = unit(-1.25, "lines"), 
    axis.text.x=element_text(angle = -90)) 


pdf(output_file, 20, 10)
plot_grid(p3, p4, ncol=2)
dev.off()









# extend the node label table to have columns that will be used for the pies
# make an additional column with the number of genes supporting the species tree bipartition
# then one column with the genes that have nothing to say
# then one column with the total number of other genes involved in conflicts
# then one column with the number of genes supporting the first top alternative
# then one column with the number of top alternatives


# see exact definitions below in the plotting section

k2$Agree  <- rep(0, length(k2$NodeP))
k2$Neutral <- rep(0, length(k2$NodeP))
k2$OtherConf <- rep(0, length(k2$NodeP))
k2$TopAlt <- rep(0, length(k2$NodeP))
k2$TopAltNum <- rep(0, length(k2$NodeP))


k3$TopAlt[1] <- 0 # no alternatives for the firt node because rooted
k3$OtherGenes[1] <- as.numeric(GenesNum) - (as.numeric(k3$tPcon[1]) + as.numeric(k3$tPconf[1]))

for (x in 2:length(k3$NodeP)) {
  if (k3$tPconf[x] == 0) {
    k3$TopAlt[x] <- 0
    k3$OtherConf[x] <- as.numeric(k3$tPconf[x]) - k3$TopAlt[x]
    k3$OtherGenes[x] <- as.numeric(GenesNum) - (as.numeric(k3$tPcon[x]) + as.numeric(k3$tPconf[x]))
  } else {
    k3$TopAlt[x] <- as.numeric(alt_dictPkey[[as.character(k3$NodeP[x])]])
    k3$OtherConf[x] <- as.numeric(k3$tPconf[x]) - k3$TopAlt[x]
    k3$OtherGenes[x] <- as.numeric(GenesNum) - (as.numeric(k3$tPcon[x]) + as.numeric(k3$tPconf[x]))
  }
}


# make a table with nodes coordinates and percentages needed for the pies

k3_per <- k3[,1:2]
k3_per$Concording <-  rep(0, length(k3_per$NodeP))
k3_per$TopConflict <-  rep(0, length(k3_per$NodeP))
k3_per$OtherConflicts <-  rep(0, length(k3_per$NodeP))
k3_per$OtherGenes <-  rep(0, length(k3_per$NodeP))

for (x in 1:length(k3$NodeP)) {
  k3_per$Concording[x] <- (as.numeric(k3$tPcon[x])*100)/GenesNum
  k3_per$TopConflict[x] <- (as.numeric(k3$TopAlt[x])*100)/GenesNum
  k3_per$OtherConflicts[x] <- (as.numeric(k3$OtherConf[x])*100)/GenesNum
  k3_per$OtherGenes[x] <- (as.numeric(k3$OtherGenes[x])*100)/GenesNum
}

# copy the table to pies_LS

pies_LS <- k3_per




#############################
# plot the pies on the tree
# same legend as in M. Johnson's script to avoid confusions and allow cross-validation
# blue= number of genes with a concording bipartition (number from the tP_con tree)
# green = number of genes with a the top alternative bipartition (maximal number from the alts file). 
# If there are more than one top alternative, count the number of genes representing one of them (chosen randomly), 
# whereas the genes representing other top alternatives are counted in the "other conflicts" part (see below)
# red = number of genes with other conflicts than the top alternative (number from the tP_con tree - maximal number from the alts file) 
# grey = number of genes that have nothing to say about that bipartition because of missing data (gene number - the three above categories)



# plot tree with node numbered like in phyparts or following R numbering

t$node_numbPT <- rep("", length(t$node_numbT))
for (x in 1:length(t$node_numbT)) {
  t$node_numbPT[x] <- paste(t$node_numbP[x], t$node_numbT[x], sep="_")
}

pPT <- ggtree(t, ladderize=T, branch.length = "none") + 
  geom_tiplab(size=2.5) + 
  xlim_tree(23) + 
  ggtitle(Title_pies) + 
  geom_nodelab(label=t$node_numbPT, size=1, vjust=0.5, hjust= -0.1)

# make pies

pies_LS_img <- nodepie(pies_LS, cols=3:6, color=c(Concording='blue', TopConflict='green', OtherConflicts='red', OtherGenes='gray'))
pies_HS_img <- nodepie(pies_HS, cols=3:6, color=c(Concording='blue', TopConflict='green', OtherConflicts='red', OtherGenes='gray'))

#plot pies

p2 <- inset(pPT, pies_LS_img, width=1.5, height=1.5, vjust=-0.3, hjust=0.5)
p3 <- inset(p2, pies_HS_img, width=1.5, height=1.5, vjust=0.5, hjust=0.5)

# wait a bit before running this
# p3







