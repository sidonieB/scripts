##############################################################################################################################
# Remove the sequences in "INTRUS" from the alignment and check if a given sequence, for instance the outgroup (OG) is there
# Can also chose to keep only the sequences in INTRUS
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
###############################################################################################################################

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

INTRUS = ["Basselinia-deplanchei-NE05_3233_S27_L001" , "Basselinia-eriostachys-NE03_2495_S16_L001" , "Basselinia-favieri-NC10_3234_S40_L001" , 
"Basselinia-gracilis-NC04_3236_S22_L001" , "Basselinia-humboldtiana-ND01_3237_S4_L001" , "Basselinia-iterata-NE11_3238_S44_L001" , 
"Basselinia-moorei-NF11_3244_S45_L001" , "Basselinia-pancheri-ND03_2494_S15_L001" , "Basselinia-porphyrea-NF03_2496_S17_L001" , 
"Basselinia-sordida-NB10_3239_S39_L001" , "Basselinia-tomentosa-NF01_3240_S5_L001" , "Basselinia-velutina-ND02_3241_S9_L001" , 
"Basselinia-vestita-NE02_3245_S10_L001" , "Cyphophoenix-elegans-SBL57_3281_S42_L001"]

x = 0

handle = open(input)
for seq_record in SeqIO.parse(handle, "fasta"):
    ID = seq_record.id
    print ID
    #ID = re.split(r"-g\d", ID)[0]
    print ID
    SEQ = str(seq_record.seq)
    SEQ = SEQ.replace("-", "")
    SEQ_d = Seq(SEQ, generic_dna)
    #print SEQ_d
    seq_record2 = SeqRecord(seq=SEQ_d, id=ID, description="")
    if "Cyphophoenix-elegans-SBL57_3281_S42_L001" in ID:        # set OG
        x = x + 1
    print x
handle.close()


if x > 0:
    print "filtering..."
    handle = open(input)
    for seq_record in SeqIO.parse(handle, "fasta"):
        ID = seq_record.id
        #ID = re.split(r"-g\d", ID)[0]
        SEQ = str(seq_record.seq)
        SEQ = SEQ.replace("-", "")
        SEQ_d = Seq(SEQ, generic_dna)
        #print SEQ_d
        seq_record2 = SeqRecord(seq=SEQ_d, id=ID, description="")
        if ID in INTRUS:            # reverse conditions if wants to filter out or keep the content of INTRUS
            with open(output, 'at') as fo:
                SeqIO.write(seq_record2, fo, "fasta")
        else:
            print "intrus en vue"
else:
    print "NO OUTGROUP!"


handle.close()
