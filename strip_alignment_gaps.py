###################################################################################################
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
# Remove columns from alignment if they have more than "threshold" percent of gaps
# Can loop on many files
####################################################################################################



#!/usr/bin/python
# -*- coding: utf-8 -*-

import re
import sys
import subprocess
from string import *
from Bio import SeqIO
import getopt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import operator
from operator import itemgetter
from Bio import AlignIO

### command: python strip_alignment_gaps.py infile.fasta outfile.fasta threshold

input = sys.argv[1]
output = sys.argv[2]
threshold = int(sys.argv[3])

alignment = AlignIO.read(open(input), "fasta")

Al_len = alignment.get_alignment_length()
tax_num = 0
for record in alignment:
   tax_num = tax_num + 1

### keep only columns whith less than threshold gaps

GOOD_COL = []
for c in range(0,Al_len,1):
   gap_num = alignment[:,c].count("-")
   gap_per = (gap_num*100)/tax_num
   if gap_per < threshold:
      GOOD_COL.append(c)

target = alignment[:,:0]
for c in GOOD_COL[0:]:
    target = target + alignment[:,c:c+1]

### write alignments in output files

AlignIO.write(target, output, "fasta")


