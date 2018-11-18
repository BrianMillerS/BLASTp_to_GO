#!/usr/bin/python

"""
Takes in a fa and alters the headerlines ">" and writes to a file.  This takes the parent= section and makes it the headerlines. 
This is done so that the drosophila header returned from BLASTp is the gene identifier and not the protein identifier, 
the GO anntoations are gene based not protein based.
"""

from Bio import SeqIO

proteome = "dmel-all-translation-r6.23.fasta"
loaded_proteome = SeqIO.parse(proteome, "fasta")

output = []
for header_line in loaded_proteome:
    full_headerline = header_line.description
    headerline_list = full_headerline.split(" ")

    # split the header and get the section with the gene name
    for split in headerline_list:
        if "parent=" in split:
            new_sequence_header = split[7:18]

    sequence = str(header_line.seq)
    output.append(">" + new_sequence_header + "\n")
    output.append(sequence + "\n")


filename = 'dmel-all-translation-r6.23_renamed_headerlines.fasta'
with open(filename, 'w+') as file:
    file.writelines(i for i in output)
