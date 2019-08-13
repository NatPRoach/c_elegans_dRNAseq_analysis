#!/usr/bin/env python2

infile = open("WS265/cds.bed")

outfile1 = open("WS265/stop_codons.bed",'w')
outfile2 = open("WS265/start_codons.bed",'w')
for line in infile:
    fields = line.strip().split()
    if fields[5] == '+':
        fields[1] = str(int(fields[2]) - 3)
    elif fields[5] == '-':
        fields[2] = str(int(fields[1]) + 3)
    outfile1.write('\t'.join(fields) + '\n')
    
    fields = line.strip().split()
    if fields[5] == '+':
        fields[2] = str(int(fields[1]) + 3)
    elif fields[5] == '-':
        fields[1] = str(int(fields[2]) - 3)
    outfile2.write('\t'.join(fields) + '\n')