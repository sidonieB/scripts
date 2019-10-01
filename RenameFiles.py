###################################################################################################
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
# Rename files based on user defined dictionnary of old:new names + modifications
####################################################################################################

import re
import sys
import os
from string import *
import getopt
import operator
from operator import itemgetter
import json


inputdir = sys.argv[1] # directory with the files to rename
output = sys.argv[2]

names_dic = {} # user-defined dictionary of old:new names


for item in os.listdir(inputdir):

   oldname = str(item)											#Dypsis-jurassic-park_S23_L001_R1_001.fastq
   paftol_id = oldname.split("_S")[0]							#Dypsis-jurassic-park
   S = oldname.split("_L001")[0]								#Dypsis-jurassic-park_S23
#   S2 = S.replace("_", "-")									#Dypsis-jurassic-park-S23
   S2 = S.split("_")[1]											#S23
   suffix = oldname.split("_L001")[1]							#_R1_001.fastq
   taxon = names_dic[paftol_id]									#Dypsis-jurassic-park-170512-3
   taxon2 = taxon.replace("_", "-")								#Dypsis-jurassic-park-170512-3
   new_name = taxon2 + "-" + S2 + suffix						#Dypsis-jurassic-park-170512-3-S23_R1_001.fastq
   rename_cmd = "mv " + oldname + " " + new_name + "\n"
   with open(output, "at") as fo:
      fo.write(rename_cmd)

