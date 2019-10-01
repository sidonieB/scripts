###################################################################################################
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
# Remove all gap of alignment: becomes a list of sequences to align again
# Useful if filtered alignment and want to realin
# Can loop on many files
####################################################################################################


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
from Bio.Alphabet import generic_dna, generic_protein


input = sys.argv[1]
output = sys.argv[2]

handle = open(input)
for seq_record in SeqIO.parse(handle, "fasta"):
    ID = seq_record.id
    SEQ = str(seq_record.seq)
    SEQ = SEQ.replace("-", "")
    SEQ_d = Seq(SEQ, generic_dna)
    #print SEQ_d
    seq_record2 = SeqRecord(seq=SEQ_d, id=ID, description="")
    with open(output, 'at') as fo:
        SeqIO.write(seq_record2, fo, "fasta")
handle.close()
