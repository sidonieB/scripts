# -*- coding: utf-8 -*-

# match name + attribute of table 1 with name + attribute of table 2
# with dictionnary of names if necessary 
# control what does not match and what is not found

import sys
import os                                       # Import necessary python modules
import re
import Bio
from collections import OrderedDict
from string import *
from Bio import SeqIO
import getopt
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from operator import itemgetter

infile1 = sys.argv[1]						# table to make a dictionary of names (e.g. published\taccepted) published_names_vs_final_species_names_used_noAmbigDup.txt
infile2 = sys.argv[2]						# table 2 (e.g. species\tGenome size) status_table_all_species.txt
infile3 = sys.argv[3]						# table 1 (e.g. species\tclimate) palm_occurence_all_last_from_Yijing.csv
outfile1 = sys.argv[4]						# table 3 (e.g. species\taccepted\tGenome size\tclimate) palm_occurence_all_last_from_Yijing_withConsStatus.csv


# Make a names dictionary

names_dict = {}

with open(infile1) as f1:
   lines1 = f1.readlines()
   for l1 in lines1[1:]:
      pub_name = l1.split("\n")[0].split("\t")[0]
      acc_name = l1.split("\n")[0].split("\t")[2]
      names_dict[pub_name] = True
      names_dict[pub_name] = acc_name

# Make an attribute dictionary based on the names dictionary and table 2 (because table 1 has multiple occurrences of the same name, if not it could be also done with table 1)
# Print a warning if two rows have a same published name (should not be the case for us)

attribute_dict = {}

with open(infile2) as f2:
   lines2 = f2.readlines()
   for l2 in lines2[1:]:
      name = l2.split("\n")[0].split("\t")[2]
      attribute = l2.split("\n")[0].split("\t",3)[3]
      #print name
      if name in names_dict:
         acc_name = names_dict[name]
         #print acc_name
         if " " in acc_name:
            acc_genus = acc_name.split(" ")[0]
            acc_species = acc_name.split(" ")[1]
            res = acc_genus + "\t" + acc_species + "\t" + acc_name + "\t" + l2
         else:
            acc_genus = "NA"
            acc_species = "NA"
            res = acc_genus + "\t" + acc_species + "\t" + acc_name + "\t" + l2
         if acc_name in attribute_dict:
            print "WARNING, name multiple times in status table; name is " + str(name) + " , accepted name is " + str(acc_name) + " and previous accepted and original names were " + str(attribute_dict[acc_name])
         else:
            attribute_dict[acc_name] = True
            attribute_dict[acc_name] = res
      else:
         print "WARNING, name in status table not in name dictionary; name is " + str(name)

# match with last table and output results

header = "Accepted_name_occ" + "\t" + "lon" + "\t" + "lat" + "\t" + "coordinateUncertaintyInMeters" + "\t" + "conservation_status_assessed_occ" + "\t" + "Acc_genus" + "\t" + "Acc_species" + "\t" + "Acc_genus_species" + "\t" + "Genus" + "\t" + "Species" + "\t" + "Genus_species" + "\t" + "Conservation_status" + "\t" + "Prediction" + "\t" + "Proba_LC" + "\t" + "Proba_nonLC"
with open(outfile1, "a") as fo1:
   fo1.write(header + "\n")
			
with open(infile3) as f3:
   lines3 = f3.readlines()
   for l3 in lines3[1:]:
      name = l3.split("\n")[0].split(",")[0]
      #print name
      attribute = l3.split("\r\n")[0].split(",",1)[1].replace(",","\t")
      if name in names_dict:
         acc_name = names_dict[name]
         if " " in acc_name:
            if acc_name in attribute_dict:
               res = name + "\t" + attribute + "\t" + attribute_dict[acc_name]
               with open(outfile1, "a") as fo1:
                  fo1.write(res)
            else:
               print "WARNING, name in occurrence table not in attribute dictionary; name is " + str(name) + " and accepted name is " + str(acc_name)
         else:
            print "WARNING, name in occurrence table has NA for accepted name; name is " + str(name)
      else:
         print "WARNING, name in occurrence table not in name dictionary; name is " + str(name)


