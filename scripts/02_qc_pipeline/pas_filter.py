#!/usr/bin/env python2

import sys
import sets
import numpy
import matplotlib
import matplotlib.pyplot as plt


#./pas_filter.py --stop_predicted    ${combined_prefix}_stringent2_predicted_stop.fa     ${combined_prefix}_stringent2_predicted_stop.bedfa     ${combined_prefix}_stringent2_has_PAS.bedfa 
fasta_in  = open(sys.argv[1])
bedfa_in  = open(sys.argv[2])
bedfa_out = open(sys.argv[3],'w')
# bedfa_out1 = open(sys.argv[3],'w')
# bedfa_out2 = open(sys.argv[3],'w')
noncanonPASTable = open("../../references/utrs/mangoneAltPAS.txt")
window = 60 
canonCount = 0
noncanonCount = 0
noPASCount = 0
altPAStable = []
canon_offsets = []
noncanon_offsets = []
canon_PAS_reads = sets.Set()
alt_PAS_reads = sets.Set()
no_PAS_reads = sets.Set()

for line in noncanonPASTable:
    fields = line.strip().split('\t')
    altPAStable.append(fields[0])

last_line = None 

for line in fasta_in:
    if line[0] == '>':
        read_id = line.strip().strip(">()+-")
        continue
    seq = line.strip().upper()
    canonFlag = False
    offset = 0
    for x in reversed(range(20,window-5)):
        subseq = seq[x:x+6]
        if subseq == "AATAAA":
            offset = window * 2  - 19 - x
            canon_offsets.append(x - window) ## TODO: Make sure this math is right
            canonFlag = True
            break
    if canonFlag:
        canonCount += 1
        canon_PAS_reads.add(read_id)
        # pasAssignmentFile.write("%s\t%s\n" %(read_id,"AATAAA"))
    else:
        noncanonFlag = False
        for altPAS in altPAStable:
            for x in reversed(range(20,window-5)):
                subseq = seq[x:x+6]
                if subseq == altPAS:
                    noncanonFlag = True
                    break
            if noncanonFlag:
                break
        if noncanonFlag:
            noncanonCount += 1
            # pasAssignmentFile.write("%s\t%s\n" %(read_id,altPAS))
            alt_PAS_reads.add(read_id)
        else:
            noPASCount += 1
            no_PAS_reads.add(read_id)
            # pasAssignmentFile.write("%s\t%s\n" %(read_id,"noPAS"))



print canonCount
print 100.*float(canonCount) / float(canonCount+noncanonCount+noPASCount)
print noncanonCount
print 100.*float(noncanonCount) / float(canonCount+noncanonCount+noPASCount)
print noPASCount
print 100.*float(noPASCount) / float(canonCount+noncanonCount+noPASCount)


for line in bedfa_in:
    fields = line.strip().split()
    if fields[3] in canon_PAS_reads or fields[3] in alt_PAS_reads:
        bedfa_out.write(line)

# fastaFile = open(sys.argv[1],'r')
# noncanonPASTable = open(sys.argv[2],'r')
# pasAssignmentFile = open(sys.argv[3],'w')
# outfiles_prefix = sys.argv[4]
#
# plotAllPAS(fastaFile,noncanonPASTable,pasAssignmentFile,outfiles_prefix)
