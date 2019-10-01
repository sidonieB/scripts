###################################################################################################
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
# Get some stats from sam depending if position on or off target (targets are defined by the input bed)
####################################################################################################

import os
import sys
import re


inputBed = sys.argv[1]
inputCov = sys.argv[2]
inputSam = sys.argv[3]
outputCov = sys.argv[4]
outputOv = sys.argv[5]

# get the on-target pos in a dict

handle = open(inputBed, "r")
lines = handle.readlines()
dictON = {}
for l in lines:
   contig = l.split("\t")[0]
   if dictON[contig] is True:
      for i in range(l.split("\t")[1], l.split("\t")[2]):
         dict[contig] = dict[contig] + "_" + "-" + i + "-"
   else:
      dict[contig] = true
      dict[contig] = str(l.split("\t")[1])
      for i in range(l.split("\t")[1], l.split("\t")[2]):
         dict[contig] = dict[contig] + "_" + "-" + i + "-"

#iterate over the coverage file and add if position is ON or OFF target

handle = open(inputCov, "r")
lines = handle.readlines()
for l in lines:
   contig = l.split("\t")[0]
   pos = l.split("\t")[1]
   pos2 = "-" + str(pos) + "-"
   cov = l.split("\t")[2]
   ONtarget = dict[contig].split("_")
   if pos in ONtarget:
      with open(outputCov, "a") as fo:
         fo.write(contig + "\t" + pos + "\t" + cov + "\tON\n")
   else:
      with open(outputCov, "a") as fo:
         fo.write(contig + "\t" + pos + "\t" + cov + "\tOFF\n")

# get overlap data from sam

with open(inputSam) as infile:
   for l in infile:
      if l.startswith("@"):
         print "header"
      else:
         read = l.split("\t")[0]
         contig = l.split("\t")[2]
         al_start = l.split("\t")[3]
         al_length = len(l.split("\t")[9])
         al_end = al_start + al_length
         ONtarget = dict[contig].split("_")
         target_overlap = 0
         for x in range(al_start, al_end):
            x = "-" + x + "-"
            if x in ONtarget:
               target_overlap = target_overlap + 1
         with open(outputOv, "a") as fo:
            fo.write(read + "\t" + contig + "\t" + target_overlap + "\n")
