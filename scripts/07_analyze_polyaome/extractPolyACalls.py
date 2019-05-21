#!/usr/bin/env python2

infile = open("../../references/polya/GSE104501_N2_TailLength.csv")
outfile = open("../../references/polya/polyALengths.txt",'w')

polyA_lengths = []

for i,line in enumerate(infile):
    if i  == 0:
        continue
    fields = line.strip().split(',')
    for x in range(10,len(fields)):
        if x == 10:
            assert fields[x][:2] == "\"["
        elif x == len(fields) - 1:
            assert fields[x][-2:] == "]\""
        polyA_lengths.append(int(fields[x].strip('\"[]')))

for polyA in polyA_lengths:
    outfile.write("%d\n"%(polyA))