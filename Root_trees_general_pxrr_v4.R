###################################################################################
# This script generates commands to run pxrr rooting function of the package phyx
# Version 4 
# 05 January 2019
# Originally written by Sidonie BELLOT
# Use and modify as you wish, but please don't hesitate to give feedback!
# on a set of trees, with one command per tree.
# It accomodates different outgroups per tree, and deals with non-monophyly of the
# user-given outgroup.
# If none of the provided outgroups is found, it will put the first provided outgroup 
# in the pxrr command. I think it is ok because it will give an error anyway with pxrr 
# (easy to check by running "grep 'not found' *PxrrRooted.tre" after you run all pxrr commands),
# but it could probably be dealt with by giving a warning here.
# The script could be easily changed to perform the rooting in R if needed.
###################################################################################

library(ape)
library(phytools)
library(gtools)


###############################################
# A few things to set before you run the script
###############################################


inputDir <- "/Users/sidob/Desktop/Palms/PAFTOL/Basselinia/Trees/gene_trees/"     # directory with the trees and the outgroup (OG) list only
inputOG <- "outgroups.txt"                                        # outgroups list (see below about how to make it)
output <- "A_pxrr_commands"                                     # name of the output (list of pxrr commands)
pattern_trees <- "_tree"                                              # whatever pattern is common to all your tree files AND ONLY TO THEM
path_to_pxrr <- "/home/rodrigo/software/phyx/src/pxrr"            # absolute path to pxrr on the computer where you will run the pxrr command

# How to make the outgroups list:
# List of the outgroups in order of priority, and with OG in clades on the same line and separated by " ; "
# The OG list has to be in quotes as follows:
# "OG1"
# "OG2"
# "OG3 ; OG4"
# etc
# If you have a lot of missing taxa, the list may contain many possible OGs
# It is ok to have only one OG in the OG list if you know that it is present in all trees




####################################################################
# You should not have to change anything below this line
# Except if you want a cleaner output file name in the pxrr command
# See comment at the very end of the script
####################################################################

setwd(inputDir)
con.trees <- list.files(".", pattern = pattern_trees, full.names = T) 
OGlist <- read.table(inputOG, stringsAsFactors = F) 
con.trees # check the input trees

OGvec <- OGlist$V1                     
OGvec     # check the input OG list                             


for (z in 1:length(con.trees)) {                    # for each tree
  
  tr <- read.tree(con.trees[[z]])
  #tr <- read.tree("RAxML_bipartitionsBranchLabels.ML.g4471.FNA.aligned_renamed.fasta.fas")
  #plot(tr)
  x <- 0
  

  for (y in 1:length(OGvec)) {          # for each OG(s), starting with the first line
        
    OGs <- strsplit(OGvec[y], " ; ")[[1]]
    
    if (length(OGs) < 2) {                          # if single OG
      if (OGs[[1]] %in% tr$tip.label) {               # check if in the tree
                                                      # if so, root the tree and stop the loop and start next tree (if not in tree, look for next OG)
        if (x < 1) {  
           final_OG <- OGs[[1]]
           x <- x + 1}
      }
    } 
    else {                                          # if multiple OG
      present_OGs <- c()
      for (o in 1:length(OGs)) {
        if (OGs[[o]] %in% tr$tip.label) {           # lists the OG that are actually present in the tree
          present_OGs <- c(present_OGs, OGs[[o]])
        } 
      }
      if (length(present_OGs) > 0) {
        if (is.monophyletic(tr, present_OGs)) {     # if those OGs are monophyletic, use them all to root the tree and stop the loop and start new tree
          if (x < 1) {
             final_OG <- present_OGs
             x <- x + 1}
        } 
        else {                      # if the OGs are not monophyletic, looks for the largest combination of those OGs that is monophyletic and use it to root the tree
          length_final_OG <- 0
          OG_combs <- do.call("c", lapply(seq_along(present_OGs), function(i) combn(present_OGs, i, FUN = list)))
          for (c in 1:length(OG_combs)) {
            if (is.monophyletic(tr, OG_combs[[c]])) {
              if (length(OG_combs[[c]]) > length_final_OG) {
                final_OGs <- OG_combs[[c]]
                length_final_OG <- length(OG_combs[[c]])
              }
            }
          }

          if (x < 1) {          
             final_OG <- final_OGs
             print(paste("WARNING OUTGROUP", present_OGs,"FOR TREE", con.trees[[z]], "IS NOT MONOPHYLETIC, ROOTING WITH THE LARGEST MONOPHYLETIC SUBGROUP AVAILABLE" ))
             plot(tr, main = con.trees[[z]])
             x <- x + 1}
        } 
     }
    }
  }
  

  final_OG_names <- final_OG[1]
  if (length(final_OG) > 1) {
  for (f in 2:length(final_OG)) {final_OG_names <- paste(final_OG_names, "," , final_OG[f], sep="" ) }}
  command <- paste(path_to_pxrr, " -t ", con.trees[[z]], " -g ", final_OG_names, " > ", sub("_tree", "", con.trees[[z]]), "_PxrrRooted.tre", sep="") # modify to get a cleaner output file name in the pxrr command if needed
  write(command, file=output, append=TRUE)
}


