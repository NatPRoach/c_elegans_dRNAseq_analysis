#!/usr/bin/env python2

import sys
import sets

inbedfile1 = open(sys.argv[1])
inbedfile2 = open(sys.argv[2])
fulloutfile = open(sys.argv[3],'w')
shortoutfile = open(sys.argv[4],'w')

full_length_reads = sets.Set()

for line in inbedfile1:
    fields = line.strip().split()
    full_length_reads.add(fields[3]) # add the read id to the list of full length reads
    
full_length = 0.0
short = 0.0
for line in inbedfile2:
    fields = line.strip().split()
    if fields[3] in full_length_reads:
        # fulloutfile.write(line)
        fulloutfile.write('\t'.join(fields[:12]) + '\n')
        full_length += 1.0
    else:
        # shortoutfile.write(line)
        shortoutfile.write('\t'.join(fields[:12]) + '\n')
        short += 1.0

total = full_length + short
print sys.argv[2]
print "Total Aligned Reads:\t%d" %(int(total))
print "Total Full Length:\t%d" %(int(full_length))
print "Total Short\t%d" %(int(short))
print "%d/%d = %f %% full length" %(int(full_length),int(total),100.0*full_length/total)
print "%d/%d = %f %% short" %(int(short),int(total),100*short/total)
