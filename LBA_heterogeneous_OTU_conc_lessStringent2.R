library(TESS)
library(phytools)
library(phangorn)
library(phylobase)
library(phyclust)

# Topology A: T1(T2,S) ; Topology B: T2(T1,S) ; Topology C: S(T1,T2)

treeH <- read.tree("MTconc_allBI.new")
NumTr <- length(treeH)
NumTr

# CP regions BI
#genes_names <- c("CP_CDS_IR" , "CP_CDS_IR_NoH" , "CP_CDS_IR_NoH_P" , "CP_CDS_IR_NoL" , "CP_CDS_IR_NoL_NoH" , "CP_CDS_IR_NoL_NoH_P" , "CP_CDS_IR_NoL_P" , "CP_CDS_IR_P" , "CP_CDS_SC" , "CP_CDS_SC_1-2" , "CP_CDS_SC_1-2_NoH" , "CP_CDS_SC_1-2_NoL" , "CP_CDS_SC_1-2_NoL_NoH" , "CP_CDS_SC_NoH" , "CP_CDS_SC_NoH_P" , "CP_CDS_SC_NoL" , "CP_CDS_SC_NoL_NoH" , "CP_CDS_SC_NoL_NoH_P" , "CP_CDS_SC_NoL_P" , "CP_CDS_SC_P" , "CP_Introns_IR" , "CP_Introns_IR_NoH" , "CP_Introns_IR_NoH_P" , "CP_Introns_IR_NoL" , "CP_Introns_IR_NoL_NoH" , "CP_Introns_IR_NoL_NoH_P" , "CP_Introns_IR_NoL_P" , "CP_Introns_IR_P" , "CP_Introns_SC" , "CP_Introns_SC_NoH" , "CP_Introns_SC_NoH_P" , "CP_Introns_SC_NoL" , "CP_Introns_SC_NoL_NoH" , "CP_Introns_SC_NoL_NoH_P" , "CP_Introns_SC_NoL_P" , "CP_Introns_SC_P" , "CP_NCS_IR" , "CP_NCS_IR_NoH" , "CP_NCS_IR_NoH_P" , "CP_NCS_IR_NoL" , "CP_NCS_IR_NoL_NoH" , "CP_NCS_IR_NoL_NoH_P" , "CP_NCS_IR_NoL_P" , "CP_NCS_IR_P" , "CP_NCS_SC" , "CP_NCS_SC_NoH" , "CP_NCS_SC_NoH_P" , "CP_NCS_SC_NoL" , "CP_NCS_SC_NoL_NoH" , "CP_NCS_SC_NoL_NoH_P" , "CP_NCS_SC_NoL_P" , "CP_NCS_SC_P")

# CP regions ML
#genes_names <- c("CP_CDS_IR" , "CP_CDS_IR_NoH" , "CP_CDS_IR_NoH_P" , "CP_CDS_IR_NoL" , "CP_CDS_IR_NoL_NoH" , "CP_CDS_IR_NoL_NoH_P" , "CP_CDS_IR_NoL_P" , "CP_CDS_IR_P" , "CP_CDS_SC" , "CP_CDS_SC_1-2" , "CP_CDS_SC_1-2_NoH" , "CP_CDS_SC_1-2_NoL" , "CP_CDS_SC_1-2_NoL_NoH" , "CP_CDS_SC_NoH" , "CP_CDS_SC_NoH_P" , "CP_CDS_SC_NoL" , "CP_CDS_SC_NoL_NoH" , "CP_CDS_SC_NoL_NoH_P" , "CP_CDS_SC_NoL_P" , "CP_CDS_SC_P" , "CP_Introns_IR" , "CP_Introns_IR_NoH" , "CP_Introns_IR_NoH_P" , "CP_Introns_IR_NoL" , "CP_Introns_IR_NoL_NoH" , "CP_Introns_IR_NoL_NoH_P" , "CP_Introns_IR_NoL_P" , "CP_Introns_IR_P" , "CP_Introns_SC" , "CP_Introns_SC_NoH" , "CP_Introns_SC_NoH_P" , "CP_Introns_SC_NoL" , "CP_Introns_SC_NoL_NoH" , "CP_Introns_SC_NoL_NoH_P" , "CP_Introns_SC_NoL_P" , "CP_Introns_SC_P" , "CP_NCS_IR" , "CP_NCS_IR_NoH" , "CP_NCS_IR_NoH_P" , "CP_NCS_IR_NoL" , "CP_NCS_IR_NoL_NoH" , "CP_NCS_IR_NoL_NoH_P" , "CP_NCS_IR_NoL_P" , "CP_NCS_IR_P" , "CP_NCS_SC" , "CP_NCS_SC_NoH" , "CP_NCS_SC_NoH_P" , "CP_NCS_SC_NoL" , "CP_NCS_SC_NoL_NoH" , "CP_NCS_SC_NoL_NoH_P" , "CP_NCS_SC_NoL_P" , "CP_NCS_SC_P")

# CP names

#T1 <- c("Trichosanthes_homophylla_SB10" , "Trichosanthes_kirilowii_SB17" , "Trichosanthes_lobata_SYS456" , "Trichosanthes_tubiflora_SYS506" , "Trichosanthes_baviensis_SYS319" , "Trichosanthes_pilosa_SYS242" , "Trichosanthes_truncata_SYS368")
#T2 <- c("Trichosanthes_tricuspidata_SB18" , "Trichosanthes_wallichiana_SB14" , "Trichosanthes_guizhouensis_SYS145")
#S <- c("Linnaeosicyos_amara_SYS507" , "Cyclanthera_pedata_SB21" , "Sechium_edule_SB9")


# mito BI
#genes_names <- c("MT_CDS" , "MT_CDS_1-2" , "MT_CDS_1-2_NoH" , "MT_CDS_1-2_NoL" , "MT_CDS_1-2_NoL_NoH" , "MT_CDS_NoH" , "MT_CDS_NoH_P" , "MT_CDS_NoL" , "MT_CDS_NoL_NoH" , "MT_CDS_NoL_NoH_P" , "MT_CDS_NoL_P" , "MT_CDS_P" , "MT_Introns" , "MT_Introns_NoH" , "MT_Introns_NoH_P" , "MT_Introns_NoL" , "MT_Introns_NoL_NoH" , "MT_Introns_NoL_NoH_P" , "MT_Introns_NoL_P" , "MT_Introns_P")

# mito ML
genes_names <- c("MT_CDS" , "MT_CDS_1-2" , "MT_CDS_1-2_NoH" , "MT_CDS_1-2_NoL" , "MT_CDS_1-2_NoL_NoH" , "MT_CDS_NoH" , "MT_CDS_NoH_P" , "MT_CDS_NoL" , "MT_CDS_NoL_NoH" , "MT_CDS_NoL_NoH_P" , "MT_CDS_NoL_P" , "MT_CDS_P" , "MT_Introns" , "MT_Introns_NoH" , "MT_Introns_NoH_P" , "MT_Introns_NoL" , "MT_Introns_NoL_NoH" , "MT_Introns_NoL_NoH_P" , "MT_Introns_NoL_P" , "MT_Introns_P")

#MITO names (T1 minimized because grade or Corallocarpus in S just to see)
T1 <- c("Trichosanthes_lobata_SYS456" , "Trichosanthes_tubiflora_SYS506" , "Trichosanthes_baviensis_SYS319" , "Trichosanthes_pilosa_SYS242")
T2 <- c("Trichosanthes_tricuspidata_SB18" , "Trichosanthes_wallichiana_SB14" , "Trichosanthes_guizhouensis_SYS145")
S <- c("Linnaeosicyos_amara_SYS507" , "Cyclanthera_pedata_SB21" , "Sechium_edule_SB9")
#S <- c("Linnaeosicyos_amara_SYS507" , "Cyclanthera_pedata_SB21" , "Sechium_edule_SB9" , "Corallocarpus_boehmii_SB19")

# nucl BI
#genes_names <- c("NC_CDS" , "NC_CDS_NoH" , "NC_CDS_NoH_P" , "NC_CDS_NoL" , "NC_CDS_NoL_NoH" , "NC_CDS_NoL_NoH_P" , "NC_CDS_NoL_P" , "NC_CDS_P" , "NC_ITS" , "NC_ITS_NoH" , "NC_ITS_NoH_P" , "NC_ITS_NoL" , "NC_ITS_NoL_NoH" , "NC_ITS_NoL_NoH_P" , "NC_ITS_NoL_P" , "NC_ITS_P" , "NC_Introns" , "NC_Introns_NoH" , "NC_Introns_NoH_P" , "NC_Introns_NoL" , "NC_Introns_NoL_NoH" , "NC_Introns_NoL_NoH_P" , "NC_Introns_NoL_P" , "NC_Introns_P" , "NC_RNA" , "NC_RNA_NoH" , "NC_RNA_NoH_P" , "NC_RNA_NoL" , "NC_RNA_NoL_NoH" , "NC_RNA_NoL_NoH_P2" , "NC_RNA_NoL_P" , "NC_RNA_P")

# nucl ML
#genes_names <- c("NC_CDS" , "NC_CDS_NoH" , "NC_CDS_NoH_P" , "NC_CDS_NoL" , "NC_CDS_NoL_NoH" , "NC_CDS_NoL_NoH_P" , "NC_CDS_NoL_P" , "NC_CDS_P" , "NC_ITS" , "NC_ITS_NoH" , "NC_ITS_NoH_P" , "NC_ITS_NoL" , "NC_ITS_NoL_NoH" , "NC_ITS_NoL_NoH_P" , "NC_ITS_NoL_P" , "NC_ITS_P" , "NC_Introns" , "NC_Introns_NoH" , "NC_Introns_NoH_P" , "NC_Introns_NoL" , "NC_Introns_NoL_NoH" , "NC_Introns_NoL_NoH_P" , "NC_Introns_NoL_P" , "NC_Introns_P" , "NC_RNA" , "NC_RNA_NoH" , "NC_RNA_NoH_P" , "NC_RNA_NoL" , "NC_RNA_NoL_NoH" , "NC_RNA_NoL_NoH_P" , "NC_RNA_NoL_P" , "NC_RNA_P")

#NUCLnames
#T1 <- c("Trichosanthes_homophylla_SB10" , "Trichosanthes_kirilowii_SB17" , "Trichosanthes_lobata_SYS456" , "Trichosanthes_pilosa_SYS242" , "Trichosanthes_tubiflora_SYS506" , "Trichosanthes_costata_SYS473" , "Trichosanthes_pedata_SYS487" , "Trichosanthes_baviensis_SYS318" , "Trichosanthes_baviensis_SYS319" , "Trichosanthes_truncata_SYS368")
#T2 <- c("Trichosanthes_tricuspidata_SB18" , "Trichosanthes_tricuspidata_SYS494" , "Trichosanthes_wallichiana_SB14" , "Trichosanthes_guizhouensis_SYS145")
#S <- c("Linnaeosicyos_amara_SYS507" , "Linnaeosicyos_amara_HE661304-JN560190-SYS507" , "Cyclanthera_pedata_SB21" , "Sechium_edule_SB9")



m <- matrix(data = NA, byrow=TRUE, nrow=NumTr, ncol=14, dimnames=list(c(genes_names), c("top" , "AvT1" , "MinT1" , "MaxT1", "AvT2" , "MinT2" , "MaxT2", "AvS" , "MinS" , "MaxS", "AvBL", "minBL", "maxBL", "SumBL")))

# test m <- matrix(data = NA, byrow=TRUE, nrow=NumTr, ncol=10, dimnames=list(c("accD", "atpA"), c("top" , "AvT1" , "MinT1" , "MaxT1", "AvT2" , "MinT2" , "MaxT2", "AvS" , "MinS" , "MaxS")))

for (z in 1:length(treeH)) {
  
  x <- treeH[[z]]
  
  #CP root
  #x <- root(x, "Corynocarpus_laevigata_NC014807")
  
  
  #MITO root
  if ("Gynostemma_burmania_SB16" %in% x$tip.label) { 
    x <- root(x, "Gynostemma_burmania_SB16") } else { 
        x <- root(x, "Luffa_aegyptiaca_SYS018") } 
  
  #NUCL root
  #if ("Corynocarpus_laevigata_AF206892-AF149001-AF479110" %in% x$tip.label) { 
  #  x <- root(x, "Corynocarpus_laevigata_AF206892-AF149001-AF479110") } else { 
  #    if ("Siratia_grosvenorii_PRJNA67121" %in% x$tip.label) { 
  #      x <- root(x, "Siratia_grosvenorii_PRJNA67121") } else { 
  #        x <- root(x, "Bryonia_dioica_SB8") } }
  
  
  All_BL <- vector(mode = "numeric", length = length(x$tip.label))
  
  for (i in 1:length(x$tip.label)) {
    tip <- x$tip.label[i]
    edgei <- which.edge(x, match(tip, x$tip.label))
    edgeiBL <- x$edge.length[edgei]
    All_BL[i] <- edgeiBL
    Anci <- Ancestors(x, match(tip, x$tip.label))
    ancAll_index <- match(Anci, Anci)
    ancAll_index2 <- ancAll_index - 1
    if (ancAll_index2 > 0) {
      for (j in 1:ancAll_index2) {
        edgej <- which.edge(x, Anci[j])
        edgejBL <- x$edge.length[edgej]
        All_BL[i] <- All_BL[i] + edgejBL
      }
    }
  }
  
  
  H<-nodeHeights(x)
  hbar2<-mean(H[sapply(length(x$tip)+2:x$Nnode,match, x$edge[,1]),1])          #from phytools blog: average BL, kind of useless
  
  hbar3<-sum(H[sapply(length(x$tip)+2:x$Nnode,match, x$edge[,1]),1])           #get the sum instead
  
  
  T1_2 <- vector(mode = "character", length = 0)
  T2_2 <- vector(mode = "character", length = 0)
  S2 <- vector(mode = "character", length = 0)
  
  for (i in T1) { 
    if (i %in% x$tip.label) {
      T1_2 <- append(T1_2, i)
    }
  }
  
  for (i in T2) { 
    if (i %in% x$tip.label) {
      T2_2 <- append(T2_2, i)
    }
  }
  
  for (i in S) { 
    if (i %in% x$tip.label) {
      S2 <- append(S2, i)
    }
  }
  
  
  if (length(T1_2) > 1) {
    AncT1 <- mrca.phylo(x, match(T1_2, x$tip.label))
    T1_all <- Descendants(x, AncT1, type = "tips") } else {T1_all <- match(T1_2, x$tip.label)}
  if (length(T2_2) > 1) {
    AncT2 <- mrca.phylo(x, match(T2_2, x$tip.label))
    T2_all <- Descendants(x, AncT2, type = "tips") } else {T2_all <- match(T2_2, x$tip.label)}
  if (length(S2) > 1) {
    AncS <- mrca.phylo(x, match(S2, x$tip.label))
    S_all <- Descendants(x, AncS, type = "tips") } else {S_all <- match(S2, x$tip.label)}
  
  T1T2 <- append(T1_2, T2_2)
  T1S <- append(T1_2, S2)
  T1T2S <- append(T1T2, S2)
  
  if (length(T1_all[[1]]) == length(T1_2)) {
    if (length(T2_all[[1]]) == length(T2_2)) {
      if (length(S_all[[1]]) == length(S2)) {
        
        AncT1T2 <- mrca.phylo(x, match(T1T2, x$tip.label))
        T1T2_all <- Descendants(x, AncT1T2, type = "tips")
        if (match(S2[1], x$tip.label) %in% T1T2_all[[1]]) {
          
          AncT1S <- mrca.phylo(x, match(T1S, x$tip.label))
          T1S_all <- Descendants(x, AncT1S, type = "tips")
          if (match(T2_2[1], x$tip) %in% T1S_all[[1]])  {
            
            m[z,1] <- "A"
            
          } else { m[z,1] <- "B" }
          
        } else { m[z,1] <- "C" }
        
        
        ancAll <- mrca.phylo(x, match(T1T2S, x$tip.label))
        T1_BL <- vector(mode = "numeric", length = length(T1))
        T2_BL <- vector(mode = "numeric", length = length(T2))
        S_BL <- vector(mode = "numeric", length = length(S))
        
        for (i in 1:length(T1_2)) {
          tip <- T1_2[i]
          edgei <- which.edge(x, match(tip, x$tip.label))
          edgeiBL <- x$edge.length[edgei]
          T1_BL[i] <- edgeiBL
          Anci <- Ancestors(x, match(tip, x$tip.label))
          ancAll_index <- match(ancAll, Anci)
          ancAll_index2 <- ancAll_index - 1
          if (ancAll_index2 > 0) {
            for (j in 1:ancAll_index2) {
              edgej <- which.edge(x, Anci[j])
              edgejBL <- x$edge.length[edgej]
              T1_BL[i] <- T1_BL[i] + edgejBL
            }
          }
        }
        
        for (i in 1:length(T2_2)) {
          tip <- T2_2[i]
          edgei <- which.edge(x, match(tip, x$tip.label))
          edgeiBL <- x$edge.length[edgei]
          T2_BL[i] <- edgeiBL
          Anci <- Ancestors(x, match(tip, x$tip.label))
          ancAll_index <- match(ancAll, Anci)
          ancAll_index2 <- ancAll_index - 1
          if (ancAll_index2 > 0) {
            for (j in 1:ancAll_index2) {
              edgej <- which.edge(x, Anci[j])
              edgejBL <- x$edge.length[edgej]
              T2_BL[i] <- T2_BL[i] + edgejBL
            }
          }
        }
        
        for (i in 1:length(S2)) {
          tip <- S2[i]
          edgei <- which.edge(x, match(tip, x$tip.label))
          edgeiBL <- x$edge.length[edgei]
          S_BL[i] <- edgeiBL
          Anci <- Ancestors(x, match(tip, x$tip.label))
          ancAll_index <- match(ancAll, Anci)
          ancAll_index2 <- ancAll_index - 1
          if (ancAll_index2 > 0) {
            for (j in 1:ancAll_index2) {
              edgej <- which.edge(x, Anci[j])
              edgejBL <- x$edge.length[edgej]
              S_BL[i] <- S_BL[i] + edgejBL
            }
          }
        }
        

        
        m[z,2] <- mean(T1_BL)
        m[z,3] <- min(T1_BL)
        m[z,4] <- max(T1_BL)
        m[z,5] <- mean(T2_BL)
        m[z,6] <- min(T2_BL)
        m[z,7] <- max(T2_BL)            
        m[z,8] <- mean(S_BL)
        m[z,9] <- min(S_BL)
        m[z,10] <- max(S_BL)
        m[z,11] <- mean(All_BL)
        m[z,12] <- min(All_BL)
        m[z,13] <- max(All_BL)
        m[z,14] <- hbar3
        
        
      } else { 
        m[z,1] <- "unresolved"
        m[z,11] <- mean(All_BL)
        m[z,12] <- min(All_BL)
        m[z,13] <- max(All_BL)
        m[z,14] <- hbar3 }
    } else {         
      m[z,1] <- "unresolved"
      m[z,11] <- mean(All_BL)
      m[z,12] <- min(All_BL)
      m[z,13] <- max(All_BL)
      m[z,14] <- hbar3 }
  } else { 
    m[z,1] <- "unresolved"
    m[z,11] <- mean(All_BL)
    m[z,12] <- min(All_BL)
    m[z,13] <- max(All_BL)
    m[z,14] <- hbar3 }
}    

write.table(m, file="MTconc_allML_BL_LessStringent.txt", sep = "\t")

