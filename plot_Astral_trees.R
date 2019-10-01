############################################################################################
# This script is to plot trees from ASTRAL, with piecharts representing Q scores or local PP
# It can also plot the actual score values instead of piecharts
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
#############################################################################################
library(ape)
library(phytools)
library(gtools)
library(ggtree)
library(cowplot)
library(ggimage)
# you probably don't need all the above packages but you certainly need ggimage, which requires an up-to-date R

setwd("../Desktop/tmp/August2018/Final_comparisons/") # set the working directory to the folder containing your tree


########### Example plot pie charts with Q scores on ASTRAL topologies



# import annotated ASTRAL tree (I like to use a tree rooted with pxrr directly after getting it from ASTRAL, before any editing)
# Then the tree has to be formated in nexus (for instance open it in Figtree and export it as nexus including annotations - don't do any other change in figtree, just open and export)
# Then it has to be edited this way (for example using Notepad++):
# replace [&label="[ by [&
# replace ]"] by ]
# replace ; by , EXCEPT for the ; at the end of the tree and in other lines of the nexus format
# Then you can import the edited nexus tree using the read.beast function:

AS_all <- read.beast("A/A_SpeciesTree_BP0_annotQ_rooted2_ed.nex")

# plot and ladderize the tree, without using the ASTRAL branch lengths
p <- ggtree(AS_all@phylo, ladderize=T, branch.length = "none") + geom_tiplab(size=2.5, hjust= -0.05) + xlim_tree(70) + ggtitle("ASTRAL all regions")

# check node labels if necessary
#p + geom_text2(aes(subset=!isTip, label=node), hjust=-.3, size = 4)

# format data to make pie charts (or any plot) of the Q scores (could do it for any other data from AS_all@data)
Q1 <- as.numeric(AS_A@data$q1) * 100
Q <- as.data.frame(Q1)
Q$Q2 <- as.numeric(AS_A@data$q2) * 100
Q$Q3 <- as.numeric(AS_A@data$q3) * 100
Q$node <- AS_all@data$node

# make barplots
# bars <- nodebar(Q, cols=1:3, position='dodge', color=c(Q1='red', Q2='cyan', Q3='gray'))
# inset(p, bars, x='node', width=8, height=2)

# make pie charts
pies <- nodepie(Q, cols=1:3, color=c(Q1='red', Q2='cyan', Q3='gray')) #, alpha=.6) # change style as needed, alpha is for transparency
inset(p, pies, width=4, height=4)#, hjust=-.6)

# pies give the quartet support: percentage of quartets agreeing: 
# with the branch (red), 
# with the second alternative RS|LO (cyan), 
# and with the last alternative RO|LS (gray).

##### plot QS to see their values
AS_A <- read.beast("A/A_SpeciesTree_BP10_annotQ_rooted2_ed.nex")

Q1 <- as.numeric(AS_A@data$q1) * 100
Q <- as.data.frame(Q1)
Q$Q2 <- as.numeric(AS_A@data$q2) * 100
Q$Q3 <- as.numeric(AS_A@data$q3) * 100
Q$node <- AS_A@data$node
node <- AS_A@data$node
allQS <- as.data.frame(node)
for (x in (1:length(allQS$node))) {allQS$QS[x] <- paste(round(Q$Q1[x]), "_", round(Q$Q2[x]), "_", round(Q$Q3[x]), sep = "")}

pP <- ggtree(AS_A@phylo, ladderize=T, branch.length = "none") %<+% allQS + 
  geom_nodelab(aes(x=branch, label=QS), vjust=-0.5, size=3) +
  geom_tiplab(size=3, hjust= -0.05) + xlim_tree(70) + 
  ggtitle("ASTRAL, plastid regions")

pdf("ASTRAL_BP10_QSvalues.pdf", 15, 15)
plot_grid(pP, pN, ncol=2)
dev.off()