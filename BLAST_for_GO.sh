#!/bin/bash

##  DEPENDENCIES  ## (make sure the following commands are in your $PATH)
# makeblastdb (version 2.2.31+)
# blastp (version 2.2.31+)

##  USAGE  ##
# ./BLAST_for_GO.sh [input_subject_proteome_fa] [input_query_proteome_fa] [output_tsv_file]
#
# for example the command "./BLAST_for_GO.sh human.fa mouse.fa mouse_mouse_best_matches.tsv" would first take the human fa and make a blast 
# database, then each mouse protein would be matched to the database and the best hit (given a evalue less than 1e-05) would be returned in 
# the file mouse_mouse_best_matches.tsv

# generate a BLASTp database from the subject proteome
makeblastdb -dbtype 'prot' -input_type 'fasta' -in $1

# run BLASTp, specify the desired outputs (this will produce a csv file that will be parsed in R)
blastp -query $2 -db $1 -out $3 -evalue 1e-05 -max_target_seqs 1 -outfmt "6 qseqid sseqid qlen slen length pident mismatch gapopen evalue nident"
#qseqid     database gene
#sseqid     query gene
#qlen       query length
#slen       subject length
#length     alighnment length
#pident     percent identity for the alignment
#mismatch   number of mismatches
#gapopen    number of gaps
#evalue     eval for the alignment
#nident     number of amino acid matches
