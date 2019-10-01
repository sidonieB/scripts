###################################################################################################
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
# Get the % of nodes in a tree that have >70% bootstrap support, for many trees: return a list with name of gene + that value
###################################################################################################################################


import re
import sys
import os
from string import *
import getopt
import operator
from operator import itemgetter
import json

inputdir = sys.argv[1] # directory with the files to parse
output = sys.argv[2]

for item in os.listdir(inputdir):
   filename = str(item)
   print filename
   fullpath = str(inputdir) + "/" + filename
   gene = filename.split("RAxML_bipartitions.")[1].split("_f")[0]			# change according to gene file name
   handle = open(fullpath, "r")
   lines = handle.readlines()
   tree = lines[0]
   li = tree.split(")")
   print(li)
   all_BP = []
   absent_BP = 0
   for x in li[1:len(li)-1]:
      BP = x.split(":")[0]
      print(x)
      if not BP:
         absent_BP = absent_BP + 1
      else:
         all_BP.append(BP)
   y = 0
   for item in all_BP:
      if int(item) >= 70:
         y = y+1
   tot_BP = len(all_BP) + absent_BP
   print tot_BP
   non_absent_BP = len(all_BP)
   print non_absent_BP
   per_BP = y*100/non_absent_BP
   with open(output, "at") as fo:
      fo.write(gene + "\t" + str(per_BP) + "\n")
   handle.close()
