############################################################################################
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
#############################################################################################

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

### Command: python FilterFastaNs.py -i infile.fa [-N xx -G xx -K xx -L xx] -o outfile.fa
### For instance: python FilterFastaNs.py -i ITS.fa -N 30 -o ITS_filtered.fa where 30 means that you want to discard all sequences with more than 30% of Ns.
### The command is flexible so that you can input one or more thresholds, if you don't specify N, G, NG and/or L, the sequences are not filtered by this criterion
### You need Python 2.7 installed
### Create a folder in which you put this script and your input file
### Then open a terminal and cd to the folder where you have your input file (a fasta file) and the script
### This script can run on many files with a loop such as: for f in *.fa; do (python FilterFastaNs.py -i $f -N 30 -o ${f/.fa}_filtered.fa); done


def ParseOpt(argv = sys.argv[1:]):		### def is used to create a function, the function ParseOpt takes the options of the command you give and interpret them to give them to the next function, FilterSeqs
   try:
      opts, args = getopt.getopt(argv,"hi:N:G:K:L:o:")			### the function getopt defines which options are expected, and where to get them. Here, the options are h, i, N, G, K, L and o, and all options have to be followed by an argument (symbolized by the ":") except h (this is why it is not followed by a ":". The options and their arguments will be found in the list "argv", which is a list of sys.argv[1:], which means a list of all what is written after the script name in the command you give)
   except getopt.GetoptError:									### if the command does not fit this scheme, the function returns the following error and quit
      print 'test.py -i <inputfile> -N <Npercent_threshold> -G <Gap_percent_threshold> -K <NandGap_percent_threshold> -L <Length_threshold> -o <outputfile>'
      sys.exit(2)
   for opt, arg in opts:										### if it is fine, the function looks in the parsed options + arguments which argument corresponds to which option
      if opt == '-h':											### if it finds the -h option, it returns how to use the script and quits
         print 'test.py -i <inputfile> -N <Npercent_threshold> -G <Gap_percent_threshold> -K <NandGap_percent_threshold> -L <Length_threshold> -o <outputfile>'
         sys.exit()
      elif opt in ("-i"):										### if it finds a -i option, then the argument that follows this option will be recognized as the input
         input = arg
      elif opt in ("-N"):										### same here an argument following -N will be recognized as the threshold for Ns
         t_N = arg
      elif opt in ("-G"):										### same for the gaps
         t_G = arg
      elif opt in ("-K"):										### same for gaps and Ns together
         t_NG = arg
      elif opt in ("-L"):										### same for length threshold
         t_L = arg
      elif opt in ("-o"):										### same for output file
         output = arg
   print 'Input file is ', input								### The next lines print to the screen the parameters recognized, so you can check if the options you specified were correctly interpreted
   print 'Output file is ', output
   try:
      print "Maximum percentage of Ns allowed set to " + t_N + "%"
   except NameError:
      print "Maximum percentage of Ns allowed set to 100%"
   try:
      print "Maximum percentage of gaps allowed set to " + t_G + "%"
   except NameError:
      print "Maximum percentage of gaps allowed set to 100%"
   try:
      print "Maximum percentage of Ns and gaps set to " + t_NG + "%"
   except NameError:
      print "Maximum percentage of Ns and gaps allowed set to 100%"
   try:
      print "Minimum sequence length set to " + t_L + "nucleotides."
   except NameError:
      print "Minimum sequence length set to 0 nucleotides."

   options = {}								### {} objects are dictionaries, here we build an options dictionary with the arguments you put or, if no arguments, with default values, for each option
   options["infile"] = True
   options["infile"] = input
   options["t_N"] = True
   try:										### This looks if you specified a threshold for Ns, and if yes, it takes it as the new threshold
      options["t_N"] = t_N
   except NameError:						### If not,
      options["t_N"] = 100					### it sets the default threshold for Ns at 100%, which means that the sequence will only be discarded if it contains more than 100% Ns, which is impossible so it will not be discarded
   options["t_G"] = True					### same for gaps
   try:
      options["t_G"] = t_G
   except NameError:
      options["t_G"] = 100
   options["t_NG"] = True					### etc.
   try:
      options["t_NG"] = t_NG
   except NameError:
      options["t_NG"] = 100
   options["t_L"] = True
   try:
      options["t_L"] = t_L
   except NameError:
      options["t_L"] = 0
   options["outfile"] = True
   options["outfile"] = output
   
   return options							### The function finishes here, it returns the dictionary of options with the new values when you input arguments



def FilterSeqs():							###This is another function that does the filtering itself using the options from the function ParseOpt
   
   options = ParseOpt()						### This is why the first thing this function does is to call (to run) the function ParseOpt
   
   infile = options["infile"]				### And using the options dictionary returned by ParseOpt, it sets the variables infile, t_N etc to the appropriate value, i.e. the ones you input in the command or the default ones if you did not input anything
   t_N = options["t_N"]
   t_G = options["t_G"]
   t_NG = options["t_NG"]
   t_L = options["t_L"]
   outfile = options["outfile"]
   
   ### Then the work really starts:
   
   handle = open(infile)			#opens the input
   n = 0							# sets a few variables that will allow to count the sequences that will be filtered out or kept, in order to get a summary and to control what happened
   g = 0
   k = 0
   l = 0
   y = 0
   z = 0
   for seq_record in SeqIO.parse(infile, "fasta"):		# go sequence by sequence
      z = z + 1											# add 1 to z so that we can count how many sequences will be processed (at the end z should be the total number of sequences in your input file)
      SEQ = seq_record.seq								# gets the first sequence item: seq_record is an item comprising the sequence id + the sequence itself. Calls "SEQ" the sequence itself
      print "processing " + seq_record.id + "..."		# prints to the screen the id of the sequence being processed
      SEQ_len = len(SEQ)								# gets the length of the sequence (i.e. the number of nucleotides)
      Ns = SEQ.count("N")								# counts the number of Ns in the sequence
      Ns_per = (Ns*100)/SEQ_len							# calculate the percentage of Ns in the sequence
      Gaps = SEQ.count("-")								# counts the number of gaps in the sequence
      Gaps_per = (Gaps*100)/SEQ_len						# calculate the percentage of gaps in the sequence
      NGaps = Ns + Gaps										# counts the number of N and gaps in the sequence
      NGaps_per = (NGaps*100)/SEQ_len						# calculate the percentage of Ns and gaps in the sequence
      print seq_record.id + " is " + str(SEQ_len) + " bp long and contains " + str(Ns) + " Ns, i.e. " + str(Ns_per) + "% Ns, " + str(Gaps) + " gaps, i.e. " + str(Gaps_per) + "% gaps, " + str(NGaps) + " Ns or gaps, i.e. " + str(NGaps_per) + "% Ns or gaps."	# prints to the screen the length and number of Ns, of gaps, of Ns and gaps and the corresponding percentages, of the sequence being processed

      if SEQ_len < int(t_L):							# if the sequence is too short,
         print seq_record.id + " is too short and will be discarded from the output file"				# prints that it is too short
         l = l + 1																						# add 1 to l so that at the end we know how many sequences were discarded AT LEAST because they were too short
      else:												# if the sequence is instead of a size superior to the threshold, continue with other thresholds
         if Ns_per > int(t_N):							# if the percentage of Ns is larger than the threshold
            print seq_record.id + " has too many Ns and will be discarded from the output file"			# print to the screen that the sequence will be discarded
            n = n + 1																					# add 1 to n so that we can count how many sequences where discarded AT LEAST because they had too many Ns
         elif Gaps_per > int(t_G):						# if the percentage of Ns is acceptable checks the percentage of gaps, and do the same
            print seq_record.id + " has too many gaps and will be discarded from the output file"
            g = g + 1
         elif NGaps_per > int(t_NG):					# if the percentage of N and the percentage of G are ok, checks the percentage of N+gaps, and do the same
            print seq_record.id + " has too many Ns and gaps and will be discarded from the output file"
            k = k + 1
         else:											# else, i.e. if the percentage of Ns and/or gaps is smaller than the thresholds, and the sequence is long enough
            print seq_record.id + " is good enough and will be kept in the output file"					# print to the screen that the sequence will be kept
            with open(outfile, 'at') as fo:																# opens the output file and call this open file "fo" just temporarily (complicated to explain)
               SeqIO.write(seq_record, fo, "fasta")														# write the sequence in the open output in the format fasta
            y = y + 1																					# add 1 to y so that we can count how many sequences will be kept (at the end y should be the total number of sequences in your output, and should be equal to z - l - n - g - k)

   handle.close()				# closes the input
   summary = []					# [] objects are lists, summary is a list of the values corresponding to the numbers of sequences discarded tu due to at least too many Ns (n), gaps (g), Ns and gaps (k), a too short length (l), the number of sequences kept (y) and the total number of sequences processed (z)
   summary.append(n)
   summary.append(g)
   summary.append(k)
   summary.append(l)
   summary.append(y)
   summary.append(z)
   return summary				# In addition to writing only good sequences in the output, the function finishes in returning the summary list, which will be used in the last function PrintSummary



def PrintSummary(results, outfile):		### This function just checks how many sequences were discarded and why
										### It needs two arguments: a result argument which will be the summary list returned by the FilterSeqs function, and the name of the output file, which we get separately (see below)


   if results[5] == (results[0] + results[1] + results[2] + results[3] + results[4]):	# checks that z indeed equal to n + g + k + l + y
      print "All sequences were parsed successfully"									# and outputs a success message (but of course there could be other problems, not spotted by this particular check)
   else:
      print "Something is wrong in this #§$ß! script!"									# or a warning
   if results[5] == 1:																	# The rest just uses the values of n, g, k, l, y and z to print to the script what happens, it is a bit convoluted because I distinguished cases were you have just one sequence and thus the sentence has to be written at the singular and not plural
      if results[0] > 0:
         print "Only one sequence was processed and it was discarded because it had too many Ns"
      elif results[1] > 0:
         print "Only one sequence was processed and it was discarded because it had too many gaps"
      elif results[2] > 0:
         print "Only one sequence was processed and it was discarded because it had too many Ns or gaps"
      elif results[3] > 0:
         print "Only one sequence was processed and it was discarded because it was too short"
      elif results[4] > 0:
         print "Only one sequence was processed and it was saved in " + str(output)
   else:
      print "For a total of " + str(results[5]) + " sequences:"
      if results[0] == 1:
         print "1 sequence was discarded because it had too many Ns"
      elif results[0] > 1:
         print str(results[0]) + " sequences were discarded because they had too many Ns"
      if results[1] == 1:
         print "1 sequence was discarded because it had too many gaps"
      elif results[1] > 1:
         print str(results[1]) + " sequences were discarded because they had too many gaps"
      if results[2] == 1:
         print "1 sequence was discarded because it had too many Ns or gaps"
      elif results[2] > 1:
         print str(results[2]) + " sequences were discarded because they had too many Ns or gaps"
      if results[3] == 1:
         print "1 sequence was discarded because it was too short"
      elif results[3] > 1:
         print str(results[3]) + " sequences were discarded because they were too short"
      if results[4] == 1:
         print "1 sequence was saved in " + str(outfile)   
      else:
         print str(results[4]) + " sequences were saved in " + str(outfile)


def main():															### This is the main function: it calls the other functions:
   filtering_results = FilterSeqs()									### First it calls the FilterSeqs function, which itself calls the ParseOpt function, ad then filters the sequences. The summary list that will be returned by this function is stored in a variable called "filtering_results"
   opts, args = getopt.getopt(sys.argv[1:],"hi:N:G:K:L:o:")		### Then the "main" function gets the name of the output file, by looking which arguments corresponds to the -o flag
   for opt, arg in opts:
      if opt in ("-o"):
         output = arg
   PrintSummary(results = filtering_results, outfile = output)		### Finally it gives both the result of FilterSeqs and the output name to the last function PrintSummary, which will return the summary to the screen
   
main()																### This runs the "main" function, which runs the other ones