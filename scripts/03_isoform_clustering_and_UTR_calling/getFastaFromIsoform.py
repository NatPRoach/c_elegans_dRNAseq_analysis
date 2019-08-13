#!/usr/bin/env python2

import sys

infile = open(sys.argv[1])
outfile= open(sys.argv[2],'w')

for line in infile:
    fields = line.strip().split()
    gene_id = fields[0]
    seq = fields[4]
    read_ids = fields[13]
    fasta_id = '>' + gene_id + '_' + read_ids
    outfile.write("%s\n" %(fasta_id))
    outfile.write("%s\n" %(seq))