import sys
from string import *
from Bio import SeqIO
import getopt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter

# Command:
# python Parse_GB.py sequence.gb sequence_CDS.fasta sequence_CDSu.fasta sequence_RNA.fasta sequence_RNAu.fasta sequence_NCS.fasta

# Arguments:
# sequence.gb: the plastome from GenBank, in GenBank format
# sequence_CDS.fasta = file name for the output containing all CDS
# sequence_CDSu.fasta = file name for the output containing all CDS only once
# sequence_RNA.fasta  = file name for the output containing all ribosomal RNAs
# sequence_RNAu.fasta  = file name for the output containing all ribosomal RNAs only once
# sequence_NCS.fasta  = file name for the output containing all intergenic spacers, with insightful names such as trnQ__rps16

# the CDSu file will include introns of protein coding genes
# the NCS file will include introns of trn (and their exon 2)

input = sys.argv[1]
output_CDS = sys.argv[2]
output_CDSu = sys.argv[3]
output_RNA = sys.argv[4]
output_RNAu = sys.argv[5]
output_NCS = sys.argv[6]

NCS_IDs = []

GPrec = "start"
G_end_Prec = 0
Fnum = 0

print input

gb = SeqIO.read(open(input, "r"), "genbank")
for f in gb.features:

   Fnum = Fnum + 1
   f_seq = f.extract(gb.seq)
   #print f.qualifiers				#uncomment to debug if you have an error: it will show you the last feature processed, you may want to check that all your features are attributed to a gene (if you are processing a custom *.gb file)
   #print f.type

   if f.type == "gene":
      try:
         f_id = f.qualifiers["gene"][0]
      except KeyError:
         f_id = f.qualifiers["label"][0]
         f_id = f_id.split(" ")[0]
      if "trn" in f_id:
         continue
      elif "rrn" in f_id:
         f_id = f_id + "___A___" + str(Fnum)
         f_rec = SeqRecord(seq=f_seq, id=f_id, description = "")
         with open(output_RNA, "a") as fo:
            SeqIO.write(f_rec, fo, "fasta")
      else:
         f_id = f_id + "___A___" + str(Fnum)
         f_rec = SeqRecord(seq=f_seq, id=f_id, description = "")
         with open(output_CDS, "a") as fo:
            SeqIO.write(f_rec, fo, "fasta")
   
   elif f.type == "rRNA":
      try:
         f_id = f.qualifiers["gene"][0] + "___B___" + str(Fnum)
      except KeyError:
         f_id = f.qualifiers["label"][0]
         f_id = f_id.split(" ")[0] + "___B___" + str(Fnum)
      f_rec = SeqRecord(seq=f_seq, id=f_id, description = "")
      with open(output_RNA, "a") as fo:
         SeqIO.write(f_rec, fo, "fasta")

   elif f.type == "CDS":
      try:
         f_id = f.qualifiers["gene"][0] + "___C___" + str(Fnum)
      except KeyError:
         f_id = f.qualifiers["label"][0]
         f_id = f_id.split(" ")[0] + "___C___" + str(Fnum)
      f_rec = SeqRecord(seq=f_seq, id=f_id, description = "")
      with open(output_CDS, "a") as fo:
         SeqIO.write(f_rec, fo, "fasta")

   elif f.type == "exon":
      try:
         f_id = f.qualifiers["gene"][0] + "___D___" + str(Fnum)
      except KeyError:
         f_id = f.qualifiers["label"][0]
         f_id = f_id.split(" ")[0] + "___D___" + str(Fnum)
      f_rec = SeqRecord(seq=f_seq, id=f_id, description = "")
      with open(output_CDS, "a") as fo:
         SeqIO.write(f_rec, fo, "fasta")
      

         
handleRNA = open(output_RNA, "rU")
l = SeqIO.parse(handleRNA, "fasta")
sortedRNA = [f for f in sorted(l, key=lambda x : x.id)]
SEQ2Prec = ""
IDPrec = ""
for s in sortedRNA:
   ID = s.id
   SEQ = s.seq
   SEQrc = SEQ.reverse_complement()
   SEQ2 = str(SEQ)
   SEQ2rc = str(SEQrc)
   if ID.split("___")[0] in IDPrec:
      if SEQ2 in SEQ2Prec:
         continue
      elif SEQ2rc in SEQ2Prec:
         continue
      else:
         with open(output_RNAu, "a") as fo:
            SeqIO.write(s, fo, "fasta")
         IDPrec = ID
         SEQ2Prec = SEQ2Prec + SEQ2
   else:
      with open(output_RNAu, "a") as fo:
         SeqIO.write(s, fo, "fasta")
      IDPrec = ID
      SEQ2Prec = SEQ2



handleCDS = open(output_CDS, "rU")
l = SeqIO.parse(handleCDS, "fasta")
sortedCDS = [f for f in sorted(l, key=lambda x : x.id)]
SEQ2Prec = ""
IDPrec = ""
for s in sortedCDS:
   ID = s.id
   SEQ = s.seq
   SEQrc = SEQ.reverse_complement()
   SEQ2 = str(SEQ)
   SEQ2rc = str(SEQrc)
   if ID.split("___")[0] in IDPrec:
      if "___A___" in ID:
         if SEQ2 in SEQ2Prec:
            continue
         elif SEQ2rc in SEQ2Prec:
            continue
         else:
            with open(output_CDSu, "a") as fo:
               SeqIO.write(s, fo, "fasta")
            IDPrec = ID
            SEQ2Prec = SEQ2Prec + SEQ2
      elif "___C___" in ID:
         if SEQ2 in SEQ2Prec:
            continue
         elif SEQ2rc in SEQ2Prec:
            continue
         else:
            if SEQ2[0:9] in SEQ2Prec:
               if SEQ2[-9:len(SEQ2)] in SEQ2Prec:
                  continue
               elif SEQ2rc[-9:len(SEQ2rc)] in SEQ2Prec:
                  continue
               else:
                  with open(output_CDSu, "a") as fo:
                     SeqIO.write(s, fo, "fasta")
                  IDPrec = ID
                  SEQ2Prec = SEQ2Prec + SEQ2
            elif SEQ2rc[0:9] in SEQ2Prec:
                continue
            else:
               with open(output_CDSu, "a") as fo:
                  SeqIO.write(s, fo, "fasta")
               IDPrec = ID
               SEQ2Prec = SEQ2Prec + SEQ2
      elif "___D___" in ID:
         if SEQ2 in SEQ2Prec:
            continue
         elif SEQ2rc in SEQ2Prec:
            continue
         else:
            with open(output_CDSu, "a") as fo:
               SeqIO.write(s, fo, "fasta")
            IDPrec = ID
            SEQ2Prec = SEQ2Prec + SEQ2
   else:
      with open(output_CDSu, "a") as fo:
         SeqIO.write(s, fo, "fasta")
      IDPrec = ID
      SEQ2Prec = SEQ2





for f in gb.features:

   if f.type == "source":
      continue
   elif f.type == "STS":
      continue
   elif f.type == "repeat_region":
      continue
   elif f.type == "misc_feature":
      continue
   elif f.type == "intron":
      continue
   else:
      try:
         G = f.qualifiers["gene"][0]
      except KeyError:
         print f.qualifiers
         print f.type
         if f.type == "tRNA":
            G = f.qualifiers["product"][0]
            G = G.split(" ")[0]
            print G
         else:
            G = f.qualifiers["label"][0]
            G = G.split(" ")[0]
      NCS_id = str(GPrec+"__"+G)
      NCS_id_rc = str(G+"__"+GPrec)
      if NCS_id_rc in NCS_IDs:
         continue
      else:
         NCS_IDs.append(NCS_id)
         G_start = f.location.nofuzzy_start
         if G_start > G_end_Prec:
            NCS_start = G_end_Prec + 1
            NCS_end = G_start - 1
            if (NCS_end - NCS_start) > 99:
               NCS_seq = gb.seq[NCS_start:NCS_end]
               NCS_rec = SeqRecord(seq=NCS_seq, id=NCS_id, description = "")
               with open(output_NCS, "a") as fo:
                  SeqIO.write(NCS_rec, fo, "fasta")

      G_end_Prec = f.location.nofuzzy_end
      GPrec = G