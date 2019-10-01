###################################################################################################
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
# Parse BLAST results
# Not ready-to-use
####################################################################################################




#!/Library/Frameworks/Python.framework/Versions/2.7/bin/python
# -*- coding: utf-8 -*-

import re
import sys
import os
from string import *
import getopt
import operator
from operator import itemgetter
import json
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

input = sys.argv[1]
outG = sys.argv[2]
outF = sys.argv[3]
outN = sys.argv[4]
outP = sys.argv[5]
out = sys.argv[6]
outW = sys.argv[7]

PRIMERS = []
hit = ""
n = 0
g = 0
f = 0
p = 0

handle = open(input, "r")
lines=handle.readlines()
for l in lines:
   primer = l.split("\t")[0]
   seqid = l.split("\r\n")[0].split("\t")[0]
   sseqid = l.split("\r\n")[0].split("\t")[1]
   sstart = l.split("\r\n")[0].split("\t")[4]
   send = l.split("\r\n")[0].split("\t")[5]
   qseq = l.split("\r\n")[0].split("\t")[6]
   sseq = l.split("\r\n")[0].split("\t")[7]
   mismatch = l.split("\r\n")[0].split("\t")[13]
   gaps = l.split("\r\n")[0].split("\t")[15]
   sstrand = l.split("\r\n")[0].split("\t")[18]
   qcovhsp = l.split("\r\n")[0].split("\t")[20]
   if primer in PRIMERS:
      if n == 0:
         if p == 0 :
            if f == 0:
               if "FLANK" in sseqid:
                  hit = hit + "\t" + seqid + "\t" + sseqid + "\t" + sstart + "\t" + send + "\t" + qseq + "\t" + sseq + "\t" + mismatch + "\t" + gaps + "\t" + sstrand + "\t" + qcovhsp + "\t"
                  Qseq = Seq(qseq)
                  Sseq = Seq(sseq)
                  SseqC = Sseq.complement()
                  #print(Qseq)
                  #print(SseqC)
                  try:
                     Tm = mt.Tm_NN(Qseq, check=True, strict=True, c_seq=SseqC, shift=0, dnac1=0.1, dnac2=0, selfcomp=False, Na=100, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
                  except ValueError:
                     with open(outW, "a") as fo:
                        fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'FLANK' + "\n")
                     Tm = 50000
                     print(Tm)
                  with open(outF, "a") as fo:
                     fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
                  hit = hit + str(Tm) + "\t"
                  f = f + 1
               elif "NOTH" in sseqid:
                  hit = hit + "\t" + seqid + "\t" + sseqid + "\t" + sstart + "\t" + send + "\t" + qseq + "\t" + sseq + "\t" + mismatch + "\t" + gaps + "\t" + sstrand + "\t" + qcovhsp + "\t"
                  Qseq = Seq(qseq)
                  Sseq = Seq(sseq)
                  SseqC = Sseq.complement()
                  #print(Qseq)
                  #print(SseqC)
                  try:
                     Tm = mt.Tm_NN(Qseq, check=True, strict=True, c_seq=SseqC, shift=0, dnac1=0.1, dnac2=0, selfcomp=False, Na=100, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
                  except ValueError:
                     with open(outW, "a") as fo:
                        fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'NOTH' + "\n")
                     Tm = 50000
                     print(Tm)
                  with open(outN, "a") as fo:
                     fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
                  hit = hit + str(Tm) + "\t"
                  n = n + 1
               elif "PLASTOME" in sseqid:
                  hit = hit + "\t" + seqid + "\t" + sseqid + "\t" + sstart + "\t" + send + "\t" + qseq + "\t" + sseq + "\t" + mismatch + "\t" + gaps + "\t" + sstrand + "\t" + qcovhsp + "\t"
                  Qseq = Seq(qseq)
                  Sseq = Seq(sseq)
                  SseqC = Sseq.complement()
                  #print(Qseq)
                  #print(SseqC)
                  try:
                     Tm = mt.Tm_NN(Qseq, check=True, strict=True, c_seq=SseqC, shift=0, dnac1=0.1, dnac2=0, selfcomp=False, Na=100, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
                  except ValueError:
                     with open(outW, "a") as fo:
                        fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'PLASTOME' + "\n")
                     Tm = 50000
                     print(Tm)
                  with open(outP, "a") as fo:
                     fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
                  hit = hit + str(Tm) + "\t"
                  p = p + 1
               else:
                  if g == 0:
                     hit = hit + "\t" + seqid + "\t" + sseqid + "\t" + sstart + "\t" + send + "\t" + qseq + "\t" + sseq + "\t" + mismatch + "\t" + gaps + "\t" + sstrand + "\t" + qcovhsp + "\t"
                     Qseq = Seq(qseq)
                     Sseq = Seq(sseq)
                     SseqC = Sseq.complement()
                     #print(Qseq)
                     #print(SseqC)
                     try:
                        Tm = mt.Tm_NN(Qseq, check=True, strict=True, c_seq=SseqC, shift=0, dnac1=0.1, dnac2=0, selfcomp=False, Na=100, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
                     except ValueError:
                        with open(outW, "a") as fo:
                           fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'GENE' + "\n")
                        Tm = 50000
                        print(Tm)
                     with open(outG, "a") as fo:
                        fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
                     hit = hit + str(Tm) + "\t"
                     g = g + 1
                  else:
                     continue
            else:
               if "FLANK" in sseqid:
                  continue
               elif "NOTH" in sseqid:
                  hit = hit + "\t" + seqid + "\t" + sseqid + "\t" + sstart + "\t" + send + "\t" + qseq + "\t" + sseq + "\t" + mismatch + "\t" + gaps + "\t" + sstrand + "\t" + qcovhsp + "\t"
                  Qseq = Seq(qseq)
                  Sseq = Seq(sseq)
                  SseqC = Sseq.complement()
                  #print(Qseq)
                  #print(SseqC)
                  try:
                     Tm = mt.Tm_NN(Qseq, check=True, strict=True, c_seq=SseqC, shift=0, dnac1=0.1, dnac2=0, selfcomp=False, Na=100, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
                  except ValueError:
                     with open(outW, "a") as fo:
                        fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'NOTH' + "\n")
                     Tm = 50000
                     print(Tm)
                  with open(outN, "a") as fo:
                     fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
                  hit = hit + str(Tm) + "\t"
                  n = n + 1
               elif "PLASTOME" in sseqid:
                  hit = hit + "\t" + seqid + "\t" + sseqid + "\t" + sstart + "\t" + send + "\t" + qseq + "\t" + sseq + "\t" + mismatch + "\t" + gaps + "\t" + sstrand + "\t" + qcovhsp + "\t"
                  Qseq = Seq(qseq)
                  Sseq = Seq(sseq)
                  SseqC = Sseq.complement()
                  #print(Qseq)
                  #print(SseqC)
                  try:
                     Tm = mt.Tm_NN(Qseq, check=True, strict=True, c_seq=SseqC, shift=0, dnac1=0.1, dnac2=0, selfcomp=False, Na=100, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
                  except ValueError:
                     with open(outW, "a") as fo:
                        fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'PLASTOME' + "\n")
                     Tm = 50000
                     print(Tm)
                  with open(outP, "a") as fo:
                     fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
                  hit = hit + str(Tm) + "\t"
                  p = p + 1
               else:
                  if g == 0:
                     hit = hit + "\t" + seqid + "\t" + sseqid + "\t" + sstart + "\t" + send + "\t" + qseq + "\t" + sseq + "\t" + mismatch + "\t" + gaps + "\t" + sstrand + "\t" + qcovhsp + "\t"
                     Qseq = Seq(qseq)
                     Sseq = Seq(sseq)
                     SseqC = Sseq.complement()
                     #print(Qseq)
                     #print(SseqC)
                     try:
                        Tm = mt.Tm_NN(Qseq, check=True, strict=True, c_seq=SseqC, shift=0, dnac1=0.1, dnac2=0, selfcomp=False, Na=100, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
                     except ValueError:
                        with open(outW, "a") as fo:
                           fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'GENE' + "\n")
                        Tm = 50000
                        print(Tm)
                     with open(outG, "a") as fo:
                        fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
                     hit = hit + str(Tm) + "\t"
                     g = g + 1
                  else:
                     continue
         else:
            if "PLASTOME" in sseqid:
               continue
            elif "NOTH" in sseqid:
               if n == 0:
                  hit = hit + "\t" + seqid + "\t" + sseqid + "\t" + sstart + "\t" + send + "\t" + qseq + "\t" + sseq + "\t" + mismatch + "\t" + gaps + "\t" + sstrand + "\t" + qcovhsp + "\t"
                  Qseq = Seq(qseq)
                  Sseq = Seq(sseq)
                  SseqC = Sseq.complement()
                  #print(Qseq)
                  #print(SseqC)
                  try:
                     Tm = mt.Tm_NN(Qseq, check=True, strict=True, c_seq=SseqC, shift=0, dnac1=0.1, dnac2=0, selfcomp=False, Na=100, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
                  except ValueError:
                     with open(outW, "a") as fo:
                        fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'NOTH' + "\n")
                     Tm = 50000
                     print(Tm)
                  with open(outN, "a") as fo:
                     fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
                  hit = hit + str(Tm) + "\t"
                  n = n + 1
               else:
                  continue
            elif "FLANK" in sseqid:
               if f == 0:
                  hit = hit + "\t" + seqid + "\t" + sseqid + "\t" + sstart + "\t" + send + "\t" + qseq + "\t" + sseq + "\t" + mismatch + "\t" + gaps + "\t" + sstrand + "\t" + qcovhsp + "\t"
                  Qseq = Seq(qseq)
                  Sseq = Seq(sseq)
                  SseqC = Sseq.complement()
                  #print(Qseq)
                  #print(SseqC)
                  try:
                     Tm = mt.Tm_NN(Qseq, check=True, strict=True, c_seq=SseqC, shift=0, dnac1=0.1, dnac2=0, selfcomp=False, Na=100, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
                  except ValueError:
                     with open(outW, "a") as fo:
                        fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'FLANK' + "\n")
                     Tm = 50000
                     print(Tm)
                  with open(outF, "a") as fo:
                     fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
                  hit = hit + str(Tm) + "\t"
                  f = f + 1
               else:
                  continue
            else:
               if g == 0:
                  hit = hit + "\t" + seqid + "\t" + sseqid + "\t" + sstart + "\t" + send + "\t" + qseq + "\t" + sseq + "\t" + mismatch + "\t" + gaps + "\t" + sstrand + "\t" + qcovhsp + "\t"
                  Qseq = Seq(qseq)
                  Sseq = Seq(sseq)
                  SseqC = Sseq.complement()
                  #print(Qseq)
                  #print(SseqC)
                  try:
                     Tm = mt.Tm_NN(Qseq, check=True, strict=True, c_seq=SseqC, shift=0, dnac1=0.1, dnac2=0, selfcomp=False, Na=100, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
                  except ValueError:
                     with open(outW, "a") as fo:
                        fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'GENE' + "\n")
                     Tm = 50000
                     print(Tm)
                  with open(outG, "a") as fo:
                     fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
                  g = g + 1
                  hit = hit + str(Tm) + "\t"
               else:
                  continue
      else:
         if "NOTH" in sseqid:
            continue
         elif "PLASTOME" in sseqid:
            if p == 0:
               hit = hit + "\t" + seqid + "\t" + sseqid + "\t" + sstart + "\t" + send + "\t" + qseq + "\t" + sseq + "\t" + mismatch + "\t" + gaps + "\t" + sstrand + "\t" + qcovhsp + "\t"
               Qseq = Seq(qseq)
               Sseq = Seq(sseq)
               SseqC = Sseq.complement()
               #print(Qseq)
               #print(SseqC)
               try:
                  Tm = mt.Tm_NN(Qseq, check=True, strict=True, c_seq=SseqC, shift=0, dnac1=0.1, dnac2=0, selfcomp=False, Na=100, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
               except ValueError:
                  with open(outW, "a") as fo:
                     fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'PLASTOME' + "\n")
                  Tm = 50000
                  print(Tm)
               with open(outP, "a") as fo:
                  fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
               hit = hit + str(Tm) + "\t"
               p = p + 1
            else:
               continue
         elif "FLANK" in sseqid:
            if f == 0:
               hit = hit + "\t" + seqid + "\t" + sseqid + "\t" + sstart + "\t" + send + "\t" + qseq + "\t" + sseq + "\t" + mismatch + "\t" + gaps + "\t" + sstrand + "\t" + qcovhsp + "\t"
               Qseq = Seq(qseq)
               Sseq = Seq(sseq)
               SseqC = Sseq.complement()
               #print(Qseq)
               #print(SseqC)
               try:
                  Tm = mt.Tm_NN(Qseq, check=True, strict=True, c_seq=SseqC, shift=0, dnac1=0.1, dnac2=0, selfcomp=False, Na=100, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
               except ValueError:
                  with open(outW, "a") as fo:
                     fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'FLANK' + "\n")
                  Tm = 50000
                  print(Tm)
               with open(outF, "a") as fo:
                  fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
               hit = hit + str(Tm) + "\t"
               f = f + 1
            else:
               continue
         else:
            if g == 0:
               hit = hit + "\t" + seqid + "\t" + sseqid + "\t" + sstart + "\t" + send + "\t" + qseq + "\t" + sseq + "\t" + mismatch + "\t" + gaps + "\t" + sstrand + "\t" + qcovhsp + "\t"
               Qseq = Seq(qseq)
               Sseq = Seq(sseq)
               SseqC = Sseq.complement()
               #print(Qseq)
               #print(SseqC)
               try:
                  Tm = mt.Tm_NN(Qseq, check=True, strict=True, c_seq=SseqC, shift=0, dnac1=0.1, dnac2=0, selfcomp=False, Na=100, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
               except ValueError:
                  with open(outW, "a") as fo:
                     fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'GENE' + "\n")
                  Tm = 50000
                  print(Tm)
               with open(outG, "a") as fo:
                  fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
               g = g + 1
               hit = hit + str(Tm) + "\t"
            else:
               continue
   else:
      with open(out, "a") as fo:
         fo.write(hit + "\n")
      PRIMERS.append(primer)
      n = 0
      g = 0
      f = 0
      p = 0
      seqid = l.split("\r\n")[0].split("\t")[0]
      sseqid = l.split("\r\n")[0].split("\t")[1]
      sstart = l.split("\r\n")[0].split("\t")[4]
      send = l.split("\r\n")[0].split("\t")[5]
      qseq = l.split("\r\n")[0].split("\t")[6]
      sseq = l.split("\r\n")[0].split("\t")[7]
      mismatch = l.split("\r\n")[0].split("\t")[13]
      gaps = l.split("\r\n")[0].split("\t")[15]
      sstrand = l.split("\r\n")[0].split("\t")[18]
      qcovhsp = l.split("\r\n")[0].split("\t")[20]
      hit = primer + "\t" + seqid + "\t" + sseqid + "\t" + sstart + "\t" + send + "\t" + qseq + "\t" + sseq + "\t" + mismatch + "\t" + gaps + "\t" + sstrand + "\t" + qcovhsp + "\t"
      Qseq = Seq(qseq)
      Sseq = Seq(sseq)
      SseqC = Sseq.complement()
      #print(Qseq)
      #print(SseqC)
      try:
         Tm = mt.Tm_NN(Qseq, check=True, strict=True, c_seq=SseqC, shift=0, dnac1=0.1, dnac2=0, selfcomp=False, Na=100, K=0, Tris=0, Mg=0, dNTPs=0, saltcorr=5)
      except ValueError:
         Tm = 50000
         print(Tm)
         if "FLANK" in sseqid:
            with open(outW, "a") as fo:
               fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'FLANK' + "\n")
         elif "NOTH" in sseqid:
            with open(outW, "a") as fo:
               fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'NOTH' + "\n")
         elif "PLASTOME" in sseqid:
            with open(outW, "a") as fo:
               fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'PLASTOME' + "\n")
         else:
            with open(outW, "a") as fo:
               fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + 'GENE' + "\n")
      if "FLANK" in sseqid:
         f = f + 1
         with open(outF, "a") as fo:
            fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
      elif "NOTH" in sseqid:
         n = n + 1
         with open(outN, "a") as fo:
            fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
      elif "PLASTOME" in sseqid:
         p = p + 1
         with open(outP, "a") as fo:
            fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
      else:
         g = g + 1
         with open(outG, "a") as fo:
            fo.write(primer + "\t" + qseq + "\t" + sseq + "\t" + str(Tm) + "\n")
      hit = hit + str(Tm) + "\t"