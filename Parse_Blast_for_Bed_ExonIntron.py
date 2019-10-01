###########################################################################################################
# Script to sort blast results depending on if good hits and make bed with matched ("exons") vs non matched ("introns") regions
# Version 1 
# 1 September 2019
# Originally written by Sidonie BELLOT (s.bellot@kew.org)
# Use and modify as you wish, but please don't hesitate to give feedback!
###########################################################################################################

import sys
from string import *
from Bio import SeqIO
import getopt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter

# Run the script as such:
# python Parse_Blast_for_Bed_ExonIntron.py All_genes.fasta blast_results.txt matching_detail.bed not_matching.txt matching_full.bed id_threshold bp_threshold

# design input and output

queries = sys.argv[1] # fasta file with all the queries
blast_res = sys.argv[2] # table output of blast following -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send qcovs evalue bitscore'
outfile = sys.argv[3] # bed file with intron and exon coordinates
outfile2 = sys.argv[4] # file with genes that did not have any good match
outfile3 = sys.argv[5] # file with genes that had a good match
id_threshold = sys.argv[6] # id% threshold (eg: 80 for 80%)
bp_threshold = sys.argv[7] # length threshold (eg: 80 for 80 bp)


# make a dictionary of query lengths
len_dict = {}
ALL_Q = []

handle = open(queries)  
for seq_record in SeqIO.parse(handle, "fasta"):
   ID = seq_record.id
   ALL_Q.append(ID)
   seq = seq_record.seq
   q_len = len(seq)
   len_dict[ID] = True
   len_dict[ID] = q_len
handle.close()


# parse results

EXONS = []
SUBJECTS = []
QUERIES = []
ALL_Q_IN_RES = []

handle = open(blast_res, "r")
lines=handle.readlines()
for l in lines:
   l2=l.split("\n")[0].split("\t")
   query=l2[0]
   h_subject=l2[1]
   h_perid=float(l2[2])
   h_ALlen=float(l2[3])
   h_s = int(l2[6])
   h_e = int(l2[7])
   if h_s < h_e:
      H_R = range(h_s, h_e,1)
   else:
      H_R = range(h_e, h_s,1)
   if query in QUERIES:
      if h_perid > float(id_threshold):
         if h_ALlen > float(bp_threshold):
            if h_subject not in SUBJECTS:
               SUBJECTS.append(h_subject)
            for i in H_R:
               if i not in EXONS:
                  EXONS.append(i)

   else:
      if len(QUERIES) > 0:
         q_len = len_dict[queryPrec]
         new_name = queryPrec
         if len(SUBJECTS) > 0:
            for sub in SUBJECTS:
               new_name = new_name + "__" + sub
            s=1
            e=0
            y=1
            z=1
            EXONS.sort()
            if min(EXONS) > 1:
               with open(outfile, "a") as fo:
                  fo.write(queryPrec + "\t" + str(1) + "\t" + str(min(EXONS)-1) + "\t"  + new_name + "___intron" + str(z) + "___" + str((min(EXONS)-1)) + "\n")
               z = z+1
               e = min(EXONS)-1
               s = min(EXONS)
            for x in EXONS:
               if x > e+1 : 
                  with open(outfile, "a") as fo:
                     fo.write(queryPrec + "\t" + str(s) + "\t" + str(e) + "\t" + new_name + "___exon" + str(y) + "___" + str(e-(s-1)) + "\n")
                     fo.write(queryPrec + "\t" + str(e+1) + "\t" + str(x-1) + "\t"  + new_name + "___intron" + str(z) + "___" + str((x-1)-e) + "\n")
                  y = y+1
                  z = z+1
                  s = x
                  e = x
               else:
                  e=x
            if x < q_len:
               with open(outfile, "a") as fo:
                  fo.write(queryPrec + "\t" + str(s) + "\t" + str(e) + "\t"  + new_name + "___exon" + str(y) + "___" + str(e-(s-1)) + "\n")
                  fo.write(queryPrec + "\t" + str(x+1) + "\t" + str(q_len) + "\t"  + new_name + "___intron" + str(z) + "___" + str(q_len-x) + "\n")
            else:
               with open(outfile, "a") as fo:
                  fo.write(queryPrec + "\t" + str(s) + "\t" + str(e) + "\t"  + new_name + "___exon" + str(y) + "___" + str(e-(s-1)) + "\n")
            with open(outfile3, "a") as fo3:
               fo3.write(queryPrec + "\t1\t" + str(q_len) + "\n")
         else:
            with open(outfile2, "a") as fo2:
               fo2.write(queryPrec + "\t" + "NO GOOD MATCH"+ "\n")

      EXONS = []
      SUBJECTS = []
      QUERIES.append(query)
      queryPrec = query
      if h_perid > float(id_threshold):
         if h_ALlen > float(bp_threshold):
            if h_subject not in SUBJECTS:
               SUBJECTS.append(h_subject)
            for i in H_R:
               if i not in EXONS:
                  EXONS.append(i)
      ALL_Q_IN_RES.append(query)

if len(QUERIES) > 0:
   q_len = len_dict[queryPrec]
   new_name = queryPrec
   if len(SUBJECTS) > 0:
      for sub in SUBJECTS:
         new_name = new_name + "__" + sub
      s=1
      e=0
      y=1
      z=1
      EXONS.sort()
      if min(EXONS) > 1:
         with open(outfile, "a") as fo:
            fo.write(queryPrec + "\t" + str(1) + "\t" + str(min(EXONS)-1) + "\t"  + new_name + "___intron" + str(z) + "___" + str((min(EXONS)-1)) + "\n")
         z = z+1
         e = min(EXONS)-1
         s = min(EXONS)
      for x in EXONS:
         if x > e+1 : 
            with open(outfile, "a") as fo:
               fo.write(queryPrec + "\t" + str(s) + "\t" + str(e) + "\t" + new_name + "___exon" + str(y) + "___" + str(e-(s-1)) + "\n")
               fo.write(queryPrec + "\t" + str(e+1) + "\t" + str(x-1) + "\t"  + new_name + "___intron" + str(z) + "___" + str((x-1)-e) + "\n")
            y = y+1
            z = z+1
            s = x
            e = x
         else:
            e=x
      if x < q_len:
         with open(outfile, "a") as fo:
            fo.write(queryPrec + "\t" + str(s) + "\t" + str(e) + "\t"  + new_name + "___exon" + str(y) + "___" + str(e-(s-1)) + "\n")
            fo.write(queryPrec + "\t" + str(x+1) + "\t" + str(q_len) + "\t"  + new_name + "___intron" + str(z) + "___" + str(q_len-x) + "\n")
      else:
         with open(outfile, "a") as fo:
            fo.write(queryPrec + "\t" + str(s) + "\t" + str(e) + "\t"  + new_name + "___exon" + str(y) + "___" + str(e-(s-1)) + "\n")
      with open(outfile3, "a") as fo3:
         fo3.write(queryPrec + "\t1\t" + str(q_len) + "\n")
   else:
      with open(outfile2, "a") as fo2:
         fo2.write(queryPrec + "\t" + "NO GOOD MATCH"+ "\n")


# list genes that had no match at all (i.e. not mentionned in the blast results)
for q in ALL_Q:
   if q in ALL_Q_IN_RES:
      print "got it"
   else:
      with open(outfile2, "a") as fo2:
         fo2.write(q + "\t" + "NO MATCH"+ "\n")



