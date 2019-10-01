###################################################################################################
# Originally written by Sidonie BELLOT - RBG Kew - s.bellot@kew.org
# Use and modify as you wish, but please don't hesitate to give feedback!
# Function to get features similar to an input feature from a database of fasta sequences, using blast
# For instance to get plastid genes from many plastomes
####################################################################################################


def BlastFeat(input,Database,output,inputPtmp,dbRtmp) :

   features_handle = open(input)
   outfile=open(output, 'a')
   for seq_record in SeqIO.parse(features_handle, "fasta"):				# for each gene
      SeqP = seq_record.seq
      SeqPid = seq_record.id
      inputP=open(inputPtmp, 'w')
      SeqIO.write(seq_record, inputP, "fasta")
      inputP.close()
      rosids_handle = open(Database)			# open the Plastome contig(s) to search
      dbR=open(dbRtmp, 'w')			# transfer the content in a new temporary file called db_R_tmp.fasta that will be the database
      dbR.write(Database)
      dbR.close()
      subprocess.call(["../Program Files/NCBI/blast-2.2.29+/bin/makeblastdb", "-in=db_R_tmp.fasta", "-dbtype=nucl"])	# make the database
      blastn_cline = NcbiblastnCommandline(query=inputPtmp, db=db_R_tmp, evalue=10, outfmt='"6 sseq"', out="blastout.txt", max_target_seqs=1, best_hit_overhang=0.1, best_hit_score_edge=0.1, word_size=11, gapopen=5, gapextend=2, reward=2, penalty=-3, strand="both", dust='"20 64 1"', perc_identity=0, xdrop_ungap=20, xdrop_gap=30, xdrop_gap_final=100)	#blast
      stdout, stderr = blastn_cline()
      result_handle = open("blastout.txt", "r")
      blast_records = result_handle.readlines()
      results = len(blast_records)
      x = 0
      for item in blast_records:
         print item
  #       if x == 0:						# get only the best HSP of the best hit - if too dissimilar, will not get the whole gene
   #         li=item.split("\t")
    #        DNArosidSeq0=li[0]
     #       li2=DNArosidSeq0.split("\n")
      #      DNArosidSeq=li2[0]
       #     DNArosidSeq2=DNArosidSeq.replace('-', '')
        #    DNArosidSeq3=Seq(DNArosidSeq2)
         #   DNArosidName = (str(SeqPid)+"______"+str(SeqRid))
          #  finalsequence = SeqRecord(seq=DNArosidSeq3, id=DNArosidName, description="")
           # SeqIO.write(finalsequence, outfile, "fasta")
            #x = x+1
      rosids_handle.close()
   outfile.close()
   features_handle.close()
