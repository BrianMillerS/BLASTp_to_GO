#!/usr/bin/python

"""
Takes in a fa and alters the headerlines ">" and writes to a file. From the uniprot fa headlines the gene name was not the first item after the ">" and thus BLAST was not using the correct seciont of the header
line. This takes the GN= section and makes it the headerlines.
"""

from Bio import SeqIO

proteome = "plasmodium_peptides_from_uniprot.fasta"
loaded_proteome = SeqIO.parse(proteome, "fasta")

output = []
for header_line in loaded_proteome:
    full_headerline = header_line.description
    headerline_list = full_headerline.split(" ")

    # split the header and get the section with the gene name
    for split in headerline_list:
        if "GN=" in split:
            new_sequence_header = split[3:]

    sequence = str(header_line.seq)

    output.append(">" + new_sequence_header + "\n")
    output.append(sequence + "\n")

filename = 'plasmodium_peptides_from_uniprot_renamed_headerlines.fasta'
with open(filename, 'w+') as file:
    file.writelines(i for i in output)
