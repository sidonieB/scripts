library(ape)
library(phytools)
library(gtools)
library(ggtree)
library(cowplot)
library(ggimage)
library(treeio)
library(gridExtra)
# you probably don't need all the above packages but you certainly need ggimage, which requires an up-to-date R

setwd("Palms/PAFTOL2/") # set the working directory to the folder containing your rooted tree


########### Example plot pie charts with Q scores on ASTRAL topologies
### Could do exactly the same with posterior probabilities, just replace "q1" by "p1" etc.

# import annotated ASTRAL tree (I like to use a tree rooted with pxrr directly after getting it from ASTRAL)

AS_0all <- read.astral("CDS_alM_r_g_c0all_trees_BP10_SpeciesTree_annotQ_rooted2.tre")
AS_40all <- read.astral("CDS_alM_r_g_c40all_trees_BP10_SpeciesTree_annotQ_rooted2.tre")
AS_70all <- read.astral("CDS_alM_r_g_c70all_trees_BP10_SpeciesTree_annotQ_rooted2.tre")

AS_0Bortho <- read.astral("CDS_alM_r_g_c0Bortho_trees_BP10_SpeciesTree_annotQ_rooted2.tre")
AS_40Bortho <- read.astral("CDS_alM_r_g_c40Bortho_trees_BP10_SpeciesTree_annotQ_rooted2.tre")
AS_70Bortho <- read.astral("CDS_alM_r_g_c70Bortho_trees_BP10_SpeciesTree_annotQ_rooted2.tre")


# plot and ladderize the tree, without using the ASTRAL branch lengths

p_0all <- ggtree(AS_0all@phylo, ladderize=T, branch.length = "none") + geom_tiplab(size=3.5, hjust= -0.05) + xlim_tree(70) + ggtitle("All genes (n=1116)")
p_40all <- ggtree(AS_40all@phylo, ladderize=T, branch.length = "none") + geom_tiplab(size=3.5, hjust= -0.05) + xlim_tree(70) + ggtitle("Genes with >40% taxa (n=1059)")
p_70all <- ggtree(AS_70all@phylo, ladderize=T, branch.length = "none") + geom_tiplab(size=3.5, hjust= -0.05) + xlim_tree(70) + ggtitle("Genes with >70% taxa (n=172)")

p_0Bortho <- ggtree(AS_0Bortho@phylo, ladderize=T, branch.length = "none") + geom_tiplab(size=3.5, hjust= -0.05) + xlim_tree(70) + ggtitle("All orthologous genes (n=445)")
p_40Bortho <- ggtree(AS_40Bortho@phylo, ladderize=T, branch.length = "none") + geom_tiplab(size=3.5, hjust= -0.05) + xlim_tree(70) + ggtitle("Orthologous genes with >40% taxa (n=418)")
p_70Bortho <- ggtree(AS_70Bortho@phylo, ladderize=T, branch.length = "none") + geom_tiplab(size=3.5, hjust= -0.05) + xlim_tree(70) + ggtitle("Orthologous genes with >70% taxa (n=120)")


# export tip list (useful to compare taxa present in different trees)

write(AS_0all@phylo$tip.label, file= "CDS_alM_r_g_c0all_trees_BP10_SpeciesTree_annotQ_rooted2_tips.txt")
write(AS_40all@phylo$tip.label, file= "CDS_alM_r_g_c40all_trees_BP10_SpeciesTree_annotQ_rooted2_tips.txt")
write(AS_70all@phylo$tip.label, file= "CDS_alM_r_g_c70all_trees_BP10_SpeciesTree_annotQ_rooted2_tips.txt")
write(AS_0Bortho@phylo$tip.label, file= "CDS_alM_r_g_c0Bortho_trees_BP10_SpeciesTree_annotQ_rooted2_tips.txt")
write(AS_40Bortho@phylo$tip.label, file= "CDS_alM_r_g_c40Bortho_trees_BP10_SpeciesTree_annotQ_rooted2_tips.txt")
write(AS_70Bortho@phylo$tip.label, file= "CDS_alM_r_g_c70Bortho_trees_BP10_SpeciesTree_annotQ_rooted2_tips.txt")


# check node labels if necessary
#p + geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size = 4)

# format data to make pie charts (or any plot) of the Q scores (could do it for any other data from AS_all@data)
# I made it into a function like that it can be done easily for many trees (clearly we could also put other parts of the script in a function, but that's a start...)
QS_calc <- function(AS_A) {
  Q1 <- as.numeric(AS_A@data$q1) * 100
  Q <- as.data.frame(Q1)
  Q$Q2 <- as.numeric(AS_A@data$q2) * 100
  Q$Q3 <- as.numeric(AS_A@data$q3) * 100
  Q$node <- AS_A@data$node
  return(Q)
}

Q_0all <- QS_calc(AS_0all)
Q_40all <- QS_calc(AS_40all)
Q_70all <- QS_calc(AS_70all)
Q_0Bortho <- QS_calc(AS_0Bortho)
Q_40Bortho <- QS_calc(AS_40Bortho)
Q_70Bortho <- QS_calc(AS_70Bortho)

# make barplots instead of pies
# bars <- nodebar(Q, cols=1:3, position='dodge', color=c(Q1='red', Q2='cyan', Q3='gray'))
# inset(p, bars, x='node', width=8, height=2)

# make pie charts (takes some time if tree is big)
pies_0all <- nodepie(Q_0all, cols=1:3, color=c(Q1='red', Q2='cyan', Q3='gray')) #, alpha=.6) # change style as needed, alpha is for transparency
pies_40all <- nodepie(Q_40all, cols=1:3, color=c(Q1='red', Q2='cyan', Q3='gray'))
pies_70all <- nodepie(Q_70all, cols=1:3, color=c(Q1='red', Q2='cyan', Q3='gray'))
pies_0Bortho <- nodepie(Q_0Bortho, cols=1:3, color=c(Q1='red', Q2='cyan', Q3='gray'))
pies_40Bortho <- nodepie(Q_40Bortho, cols=1:3, color=c(Q1='red', Q2='cyan', Q3='gray'))
pies_70Bortho <- nodepie(Q_70Bortho, cols=1:3, color=c(Q1='red', Q2='cyan', Q3='gray'))

p2_0all <- inset(p_0all, pies_0all, width=0.1, height=0.1)#, hjust=-.6) # change size if pies too big/small (better check first how it looks in the exported file)
p2_40all <- inset(p_40all, pies_40all, width=0.1, height=0.1)
p2_70all <- inset(p_70all, pies_70all, width=0.1, height=0.1)
p2_0Bortho <- inset(p_0Bortho, pies_0Bortho, width=0.1, height=0.1)
p2_40Bortho <- inset(p_40Bortho, pies_40Bortho, width=0.1, height=0.1)
p2_70Bortho <- inset(p_70Bortho, pies_70Bortho, width=0.1, height=0.1)

#export tree to pdf
pdf("CDS_alM_r_g_ALLtrees_BP10_SpeciesTree_annotQ_rooted2_QSpies.pdf", 30, 60) # change size and file name as required (eg for one tree with 200 taxa: 15, 45)
# for one plot:
#plot(p2)
#for many plots in one file:
grid.arrange(p2_0all, p2_40all, p2_70all, p2_0Bortho, p2_40Bortho, p2_70Bortho, nrow=2)
dev.off()


# pies give the quartet support: percentage of quartets agreeing: 
# with the branch (red), 
# with the second alternative RS|LO (cyan), 
# and with the last alternative RO|LS (gray).
# see also DiscoVista software to figure out what are second and third alternatives respectively for a given branch

##### plot QS to see their values
AS_A <- read.astral("A_trees_BP10_SpeciesTree_annotQ_rooted2.tre") # same tree as above, no need to read again if did pies first

# make a label with the three QS separated by underscores
Q1 <- as.numeric(AS_A@data$q1) * 100
Q <- as.data.frame(Q1)
Q$Q2 <- as.numeric(AS_A@data$q2) * 100
Q$Q3 <- as.numeric(AS_A@data$q3) * 100
Q$node <- AS_A@data$node
node <- AS_A@data$node
allQS <- as.data.frame(node)
for (x in (1:length(allQS$node))) {allQS$QS[x] <- paste(round(Q$Q1[x]), "_", round(Q$Q2[x]), "_", round(Q$Q3[x]), sep = "")}

# create plot (change label size etc as necessary after trying first)
p3 <- ggtree(AS_A@phylo, ladderize=T, branch.length = "none") %<+% allQS + 
  geom_nodelab(aes(x=branch, label=QS), vjust=-0.5, size=3) +
  geom_tiplab(size=3, hjust= -0.05) + xlim_tree(70) + 
  ggtitle("ASTRAL tree with quartet scores")

# export to pdf
pdf("ASTRAL_BP10_QSvalues.pdf", 15, 15) # change size and file name as required
plot(p3)
dev.off()