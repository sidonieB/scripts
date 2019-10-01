###################################################################################################
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
# Screen trees for a list of outgroups and report the tree name + the most distant outgroup found
# Outgroup list has to be ordered from the most distant (clade of) outgroups to the closest ones to the ingroup
# Just ask me
# I think I replaced it by a new R script
################################################################################################################


import sys
from string import *
from Bio import SeqIO
import getopt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter
import os
import re


inputdir = sys.argv[1]	# tree folder
inputTax = sys.argv[2]	# OG list
output = sys.argv[3]	# OG table with OG for each tree

_nsre = re.compile('([0-9]+)')
def natural_sort_key(s):
    return [int(text) if text.isdigit() else text.lower()
            for text in re.split(_nsre, s)] 

OG_tab = []

handleTax = open(inputTax, "r")
linesT = handleTax.readlines()

f_list = os.listdir(inputdir)
for input in f_list:
   x = 0
   TAX = []
   input2 = "nex/" + str(input)
   handle = open(input2, "r")
   gene = str(input)
   lines = handle.readlines()
   for l in lines:
      if "ntax" in l:
         NTAX = l.split("ntax=")[1].split(";")[0]

   for l in lines:
      if x < (int(NTAX) + 4):
         if x > 3:
            l2 = l.split("\t")[1].split("\r\n")[0]
            TAX.append(l2)
            x = x + 1
         else:
            x = x + 1
   print TAX
   handle.close()

   y = 0

   for lT in linesT:
      if y < 1 :
         OG = lT.split("\n")[0]
         print OG
         if ";" in OG:
            OG1 = OG.split(" ; ")[0]
            OG2 = OG.split(" ; ")[1]
            if OG1 in TAX:
               if OG2 in TAX:
                  res = gene + "\t" + OG1 + ";" + OG2 + "\n"
                  OG_tab.append(res)
                  y = y + 1
               else:
                  print "keep looking"
            else:
               print "keep looking"
         else:
            if OG in TAX:
               res = gene + "\t" + OG + "\n"
               OG_tab.append(res)
               y = y + 1
            else:
               print "keep looking"
      else:
         print "found OG"
   if y < 1 :
      res = gene + "\tNOT-FOUND\n"
      OG_tab.append(res)
   else :
      print "OG was found"
handleTax.close()

OG_tab.sort(key=natural_sort_key)
for item in OG_tab:
   with open(output, "a") as fo:
      fo.write(item)