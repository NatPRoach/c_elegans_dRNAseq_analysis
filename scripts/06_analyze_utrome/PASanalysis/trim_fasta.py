#!/usr/bin/env python2
## Takes in a fasta file and returns a fasta trimmed to a specific window size
##
import sys

window = 120 
infile = open(sys.argv[1])
outfile = open(sys.argv[2],'w')
seq_dict = {}
for line in infile:
    if line[0] == '>':
        name = line.strip().strip('>')
        seq_dict[name] = ""
    else:
        seq_dict[name] += line.strip()
for name in seq_dict:
    outfile.write(">" + name + "\n")
    outfile.write("%s\n" %(seq_dict[name][-window:].upper()))