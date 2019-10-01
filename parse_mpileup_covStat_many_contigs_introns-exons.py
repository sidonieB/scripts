###################################################################################################
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
# Get some stats from samtools mpileup output
####################################################################################################


import os
import sys
import re
import numpy

# works with Python2.x, may need to be adapted for Python 3 because range returns an iterator, so would need to do list(range())
# example command:
# python parse_mpileup_covStat_many_contigs_introns-exons.py sample_BWA_parsed_sorted_pileup.txt sample_Allintrons.gff sample_ReadDepthStats.txt

input = sys.argv[1]      # mpileup all genes (*_BWA_parsed_sorted_pileup.txt)
input2 = sys.argv[2]      # intron coordinates all genes (*_Allintrons.gff)
output = sys.argv[3]     # will be created automatically


with open(output, "a") as fo:
   fo.write("Contig\tMean exon read depth\tMin exon read depth\tMax exon read depth\tExon read depth std\tMean intron read depth\tMin intron read depth\tMax intron read depth\tIntron read depth std\n")

int_dict = {}

with open(input2, "r") as f:
   for line in f:
      contig = line.split("\t")[0] + "_"
      start = int(line.split("\t")[3])
      end = int(line.split("\t")[4])
      if end < start:
         end2 = start
         start2 = end
         end = end2
         start = start2
      pos_lst = range(start, end+1)
      pos_lst2 = ["_" + str(i) + "_" for i in pos_lst]
      if int_dict.has_key(contig):
         ls_tmp = int_dict[contig]
         [ls_tmp.append(x) for x in pos_lst2]
         int_dict[contig] = ls_tmp
      else:
         int_dict[contig] = True
         int_dict[contig] = pos_lst2

#print int_dict


CONTIGS = []
COV_int = []
COV_ex = []
contigPrec = ""
INT_POS = []

with open(input, "r") as f:
   for line in f:
      contig = line.split("\t")[0]
      contig2 = contig + "_"
      contig3 = contig.split("_L001-")[1]
      cov = int(line.split("\t")[3])
      pos = "_" + str(line.split("\t")[1]) + "_"
      if contig2 in CONTIGS:
         if pos in INT_POS:
            COV_int.append(cov)
         else:
            COV_ex.append(cov)
      else:
         CONTIGS.append(contig2)
         if COV_int:
            cov_int_mean = str(numpy.mean(COV_int))
            cov_int_min = str(numpy.min(COV_int))
            cov_int_max = str(numpy.max(COV_int))
            cov_int_sd = str(numpy.std(COV_int))
         else:
            cov_int_mean = "na"
            cov_int_min = "na"
            cov_int_max = "na"
            cov_int_sd = "na"
         if COV_ex:
            cov_ex_mean = str(numpy.mean(COV_ex))
            cov_ex_min = str(numpy.min(COV_ex))
            cov_ex_max = str(numpy.max(COV_ex))
            cov_ex_sd = str(numpy.std(COV_ex))
         else:
            cov_ex_mean = "na"
            cov_ex_min = "na"
            cov_ex_max = "na"
            cov_ex_sd = "na"
         with open(output, "a") as fo:
            fo.write(contigPrec+"\t"+cov_ex_mean+"\t"+cov_ex_min+"\t"+cov_ex_max+"\t"+cov_ex_sd+"\t"+cov_int_mean+"\t"+cov_int_min+"\t"+cov_int_max+"\t"+cov_int_sd+"\n")
         contigPrec = contig3
         COV_int = []
         COV_ex = []
         if int_dict.has_key(contig2):
            INT_POS = int_dict[contig2]
         else:
            INT_POS = []
         if pos in INT_POS:
            COV_int.append(cov)
         else:
            COV_ex.append(cov)
