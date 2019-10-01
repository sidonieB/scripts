###################################################################################################
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
# Rename sequences in a fasta file according to user-defined dictionnary of old:new names
####################################################################################################


import sys
from string import *
from Bio import SeqIO
import getopt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter

# command from the terminal: 
# for f in *.fasta; do (python Fasta_Rename.py $f ${f/.fasta}_renamed.fasta); done

input = sys.argv[1]
output = sys.argv[2]


names = {} # user-defined dictionnary of old:new names



handle = open(input)  
for seq_record in SeqIO.parse(handle, "fasta"):
   ID = seq_record.id
   seq = seq_record.seq
   #des = seq_record.description
   #length = len(seq)
   ID2 = names[ID.split("_")[0]]
   #ID2 = str(ID) + "__" + str(length)
   finalsequence = SeqRecord(seq=seq, id=ID2, description="")
   with open(output, "at") as fo:
      SeqIO.write(finalsequence, fo, "fasta")

handle.close()