###########################################################################################################
# Script to classify genes according to how they support each node of the tree and according to user-defined categories
# The support classification is based on the output from phyparts (S. Smith)
# Version 3 
# 16 March 2019
# Originally written by Sidonie BELLOT (s.bellot@kew.org)
# Use and modify as you wish, but please don't hesitate to give feedback!
###########################################################################################################

import re
import sys
import os
from string import *
import getopt
import operator
from operator import itemgetter


# The script has to be run from the folder containing the phyparts outputs
# command to run from the terminal:
# python Phyparts_parsing.py gene_names_cat node_key out0 out70 outfile

# For instance:
# python Phyparts_parsing_v3.py genes_names_cat.txt phyparts_A_BP0.node.key phyparts_A_BP0 phyparts_A_BP70 RAxML_bestTree_A_conc_part_rooted2_gene_summary.txt

# gene_names_cat is a list of all gene names including their path as they are displayed in the concord or conflict phyparts output files, and of their category (for instance Plastid, Nuclear). One gene and its category per line, tab separated such as gene	category
# You can obtain it by running something like "echo gene_trees/RAxML_bipartitions.* | sed 's/ /\n/g' >> gene_names.txt" from the folder where you ran phyparts
# But you will need to input the categories yourself, in excel for instance
# node_key is the phyparts output file comprising the code according which phyparts numbered the nodes. It is the file ending with ".node.key"
# Ensure that node key is the same for the 0% and 70% output of a same species tree (using diff for instance) and just provide one of them.
# out0 is the common part (prefix) of the names of the result files from the phyparts analysis without bootstrap threshold
# out70 is the common part (prefix) of the names of the result files from the phyparts analysis with bootstrap threshold -s 70





# design input and output
All_genes = sys.argv[1] 
node_key = sys.argv[2]
out0 = sys.argv[3]
out70 = sys.argv[4]
outfile = sys.argv[5]


# Create a dictionary with the node numbers as keys and the number of genes representing the highest alternative
# if there is a tie, use the number of gene representing one of them, not all of them
# the number to get is the one in brackets in the alts file: alt for 1 (3): xxxxxx
# Do this for the analysis taking all bipartitions into account and separately for the analysis after collapsing bipartitions with BP% < 70

alt_num_dict0 = {}
alt_num_dict70 = {}
alts0 = out0 + ".log"
alts70 = out70 + ".log"
NODES = []
NUMS = [0]
PrevNode = ""

handle = open(alts0, "r")
lines=handle.readlines()
for l in lines:
   if "alt for " in l:
      num = l.split("(")[1].split(")")[0]
      node = "_" + str(l.split("alt for ")[1].split(" ")[0]) + "_"
      if node in NODES:
         NUMS.append(num)
      else:
         alt_num_dict0[PrevNode] = True
         alt_num_dict0[PrevNode] = max(NUMS)
         PrevNode = node
         NUMS = []
         NODES.append(node)
         NUMS.append(num)
alt_num_dict0[PrevNode] = True
alt_num_dict0[PrevNode] = max(NUMS)
handle.close()

NODES = []
NUMS = [0]
PrevNode = ""

handle = open(alts70, "r")
lines=handle.readlines()
for l in lines:
   if "alt for " in l:
      num = l.split("(")[1].split(")")[0]
      node = "_" + str(l.split("alt for ")[1].split(" ")[0]) + "_"
      if node in NODES:
         NUMS.append(num)
      else:
         alt_num_dict70[PrevNode] = True
         alt_num_dict70[PrevNode] = max(NUMS)
         PrevNode = node
         NUMS = []
         NODES.append(node)
         NUMS.append(num)
alt_num_dict70[PrevNode] = True
alt_num_dict70[PrevNode] = max(NUMS)
handle.close()


# print dictionaries to check
# print alt_num_dict0
# print alt_num_dict70
# print alt_num_dict70["_12_"]

# Create a dictionary with the node numbers + alternative index as keys and the names of the genes representing the alternative as values
# This dictionary contains only the alternatives that are a top alternative, if the "top" is represented only by one gene, the dictionary does not include it
# The alternative index should be 1, or higher if there are ties, in which case the first tie listed in the file gets the smaller index etc.

alt_gen_dict0 = {}
alt_gen_dict70 = {}
x = 0
NODES = []

handle = open(alts0, "r")
lines=handle.readlines()
for l in lines:
   if "alt for " in l:
      num = l.split("(")[1].split(")")[0]
      node = "_" + str(l.split("alt for ")[1].split(" ")[0]) + "_"
      if node in NODES:
         if num == alt_num_dict0[node]:
            if int(num) > 1:
               x = x + 1
               genes = l.split(" ", 5)[5].replace(" ", ",")
               key_name = node + str(x)
               alt_gen_dict0[key_name] = True
               alt_gen_dict0[key_name] = genes
            else:
               print num
      else:
         NODES.append(node)
         x = 0
         if num == alt_num_dict0[node]:
            if int(num) > 1:
               x = x + 1
               genes = l.split(" ", 5)[5].replace(" ", ",")
               key_name = node + str(x)
               alt_gen_dict0[key_name] = True
               alt_gen_dict0[key_name] = genes
            else:
               print num
handle.close()

x = 0
NODES = []

handle = open(alts70, "r")
lines=handle.readlines()
for l in lines:
   if "alt for " in l:
      num = l.split("(")[1].split(")")[0]
      node = "_" + str(l.split("alt for ")[1].split(" ")[0]) + "_"
      if node in NODES:
         if num == alt_num_dict70[node]:
            if int(num) > 1:
               x = x + 1
               genes = l.split(" ", 5)[5].replace(" ", ",")
               key_name = node + str(x)
               alt_gen_dict70[key_name] = True
               alt_gen_dict70[key_name] = genes
            else:
               print num
      else:
         NODES.append(node)
         x = 0
         if num == alt_num_dict70[node]:
            if int(num) > 1:
               x = x + 1
               genes = l.split(" ", 5)[5].replace(" ", ",")
               key_name = node + str(x)
               alt_gen_dict70[key_name] = True
               alt_gen_dict70[key_name] = genes
            else:
               print num
handle.close()



#print alt_gen_dict0
#print alt_gen_dict70

# Create a dictionary of the genes and their category

genes_cat = {}

handle = open(All_genes, "r")
lines=handle.readlines()
for l in lines:
   gene = l.split("\t")[0]
   cat = l.split("\t")[1].split("\r\n")[0]
   genes_cat[gene] = True
   genes_cat[gene] = cat
handle.close()




header = "node\tgene\tcategory\tag0\tag70\tc0\tc70\tclassified\tneutral\tstatus\n"
with open(outfile, "a") as fo:
   fo.write(header)

handle = open(node_key, "r")
lines=handle.readlines()
for l in lines[1:]:										# get each node from the key, except node 0
   node = str(l.split(" ")[0])
   print(node)
   concon70 = out70 + ".concord.node." + node			# for each node, get the corresponding files of genes conflicting and concording at 0 and 70% 
   conflict70 = out70 + ".conflict.node." + node
   concon0 = out0 + ".concord.node." + node
   conflict0 = out0 + ".conflict.node." + node
   AGREE70 = []
   CONFLICT70 = []
   AGREE0 = []
   CONFLICT0 = []
   CLASSIFIED = []
   REST = []
   try:
      handle2 = open(concon70, "r")						# open the file with the genes which agree even after their topology was collapsed under 70% support 
      lines2=handle2.readlines()
      for l2 in lines2:
         gene = l2.split("\n")[0]
         gene = gene.split("\r")[0]
         AGREE70.append(gene)							# put them in the list AGREE70
         CLASSIFIED.append(gene)
      handle2.close()
   except IOError:
      pass
   try:
      handle2 = open(conflict70, "r")						# open the file with the genes which conflict even after their topology was collapsed under 70% support
      lines2=handle2.readlines()
      for l2 in lines2:
         gene = l2.split("\n")[0]
         gene = gene.split("\r")[0]
         CONFLICT70.append(gene)							# put them in the list CONFLICT70
         CLASSIFIED.append(gene)
      handle2.close()
   except IOError:
      pass
   try:
      handle2 = open(concon0, "r")							# open the file with the genes which agree when their topology is not collapsed under 70% support, and also potentially when it is collapsed
      lines2=handle2.readlines()
      for l2 in lines2:
         gene = l2.split("\n")[0]
         gene = gene.split("\r")[0]
         AGREE0.append(gene)							# put them in the list AGREE0
         CLASSIFIED.append(gene)
      handle2.close()
   except IOError:
      pass
   try:
      handle2 = open(conflict0, "r")							# open the file with the genes which conflict when their topology is not collapsed under 70% support, and also potentially when it is collapsed
      lines2=handle2.readlines()
      for l2 in lines2:
         gene = l2.split("\n")[0]
         gene = gene.split("\r")[0]
         CONFLICT0.append(gene)							# put them in the list CONFLICT0
         CLASSIFIED.append(gene)
      handle2.close()
   except IOError:
      pass
   handle2 = open(All_genes, "r")							# open the file with all the gene names
   lines2=handle2.readlines()
   for l2 in lines2:
      gene = l2.split("\n")[0]
      gene = gene.split("\r")[0].split("\t")[0]
      gene_cat = genes_cat[gene]						# get the category of the gene from the dictionary 
      if gene in AGREE70:									# check in which lists the gene is present and code it by changing the values of the variables ag_c0 etc to 0 (absent) or 1 (present)
         gene_ag70 = str(1)
      else:
         gene_ag70 = str(0)
      if gene in AGREE0:
         gene_ag0 = str(1)
      else:
         gene_ag0 = str(0)
      if gene in CONFLICT70:
         gene_c70 = str(1)
      else:
         gene_c70 = str(0)
      if gene in CONFLICT0:
         gene_c0 = str(1)
      else:
         gene_c0 = str(0)
      if gene in CLASSIFIED:								# checks if the gene is in any list at all: if not, code gene_neutral as 1
         gene_class = str(1)
         gene_neutral = str(0)
      else:
         gene_class = str(0)
         gene_neutral = str(1)
         #print(gene)
# check the final gene codes and deduce the gene status from it:
      if "_"+gene_ag0+"_"+gene_ag70+"_"+gene_c0+"_"+gene_c70+"_"+gene_class+"_"+gene_neutral+"_" in "_0_0_0_0_0_1_":
         status = "ABSENT_BIP"
      elif "_"+gene_ag0+"_"+gene_ag70+"_"+gene_c0+"_"+gene_c70+"_"+gene_class+"_"+gene_neutral+"_" in "_0_0_1_0_1_0_":
         status = "CONFLICT_LS"
      elif "_"+gene_ag0+"_"+gene_ag70+"_"+gene_c0+"_"+gene_c70+"_"+gene_class+"_"+gene_neutral+"_" in "_0_0_0_1_1_0_":
         status = "CONFLICTonly70"
      elif "_"+gene_ag0+"_"+gene_ag70+"_"+gene_c0+"_"+gene_c70+"_"+gene_class+"_"+gene_neutral+"_" in "_0_0_1_1_1_0_":
         status = "CONFLICT_HS"
      elif "_"+gene_ag0+"_"+gene_ag70+"_"+gene_c0+"_"+gene_c70+"_"+gene_class+"_"+gene_neutral+"_" in "_1_0_0_0_1_0_":
         status = "AGREE_LS"
      elif "_"+gene_ag0+"_"+gene_ag70+"_"+gene_c0+"_"+gene_c70+"_"+gene_class+"_"+gene_neutral+"_" in "_0_1_0_0_1_0_":
         status = "AGREEonly70"
      elif "_"+gene_ag0+"_"+gene_ag70+"_"+gene_c0+"_"+gene_c70+"_"+gene_class+"_"+gene_neutral+"_" in "_1_1_0_0_1_0_":
         status = "AGREE_HS"
      else:
         status = "UNCLEAR"

# in theory there should be no gene with status UNCLEAR, CONFLICTonly70 or AGREEonly70

# Change gene status if representing a top alternative of the HS analysis and/or of the LS analysis

      if "CONFLICT_LS" in status:
         x = 0
         for item in alt_gen_dict0:
            x = x + 1
            alt_key = "_" + str(node) + "_" + str(x)
            if alt_key in alt_gen_dict0:
               if gene in alt_gen_dict0[alt_key]:
                  status = "CONFLICT_LS" + "_tLS_" + str(x)
                  break                                        # get only the first top in which the gene is represented
      elif "CONFLICT_HS" in status:
         x = 0
         for item in alt_gen_dict0:
            x = x + 1
            alt_key = "_" + str(node) + "_" + str(x)
            if alt_key in alt_gen_dict0:
               if gene in alt_gen_dict0[alt_key]:
                  status = "CONFLICT_HS" + "_tLS_" + str(x)
                  break
         x = 0
         for item in alt_gen_dict70:
            x = x + 1
            alt_key = "_" + str(node) + "_" + str(x)
            if alt_key in alt_gen_dict70:
               if gene in alt_gen_dict70[alt_key]:
                  status = status + "_tHS_" + str(x)
                  break
      else:
         print "PAS TOP"

# could modify above to get already status for pies so that don't need to do to much in R
# Export gene status and other info

      gene_res = node + "\t" + gene + "\t" + gene_cat + "\t" + gene_ag0 + "\t" + gene_ag70 + "\t" + gene_c0 + "\t" + gene_c70 + "\t" + gene_class + "\t" + gene_neutral + "\t" + status + "\n"
      with open(outfile, "a") as fo:
         fo.write(gene_res)
handle.close()
   
