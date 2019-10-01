###################################################################################################
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
# script takes elements from list 2 and print them in output if they are not in list 1
####################################################################################################



# -*- coding: utf-8 -*-



import sys
import os                                       # Import necessary python modules
import re
import Bio
from collections import OrderedDict
from string import *
from Bio import SeqIO
import getopt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter

infile1 = sys.argv[1]						# list1
infile2 = sys.argv[2]						# list2
outfile = sys.argv[3]


BAD = []

with open(infile1) as f1:
   lines1 = f1.readlines()
   for l1 in lines1:
      item = l1.split("\n")[0]
      BAD.append(item)

with open(infile2) as f2:
   lines2 = f2.readlines()
   for l2 in lines2:
      item = l2.split("\n")[0]
      if item in BAD:
         print "We don't want " + str(item)
#         with open(outfile, "a") as fo:
#            fo.write(str(item) + "\tIS PRESENT\n")
      else:
         print str(item) + " IS GOOD!"
         item2 = str(item).split("trees/rooted/RAxML_bipartitions.")[1].split("_tree")[0] + ".fasta"
         with open(outfile, "a") as fo:
            fo.write("cp " + str(item2) + " AgreeAstralBP70\n")

