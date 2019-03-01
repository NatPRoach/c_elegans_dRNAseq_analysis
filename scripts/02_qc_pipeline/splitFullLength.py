#!/usr/bin/env python2

import sys
import sets
import pysam

inbedfile = open(sys.argv[1])
inbamfile = pysam.AlignmentFile(sys.argv[2],'rb')
fulloutfile = pysam.AlignmentFile(sys.argv[3],'wb',template = inbamfile)
shortoutfile = pysam.AlignmentFile(sys.argv[4],'wb',template = inbamfile)

full_length_reads = sets.Set()

for line in inbedfile:
    fields = line.strip().split()
    full_length_reads.add(fields[3]) # add the read id to the list of full length reads
    
full_length = 0.0
short = 0.0
for read in inbamfile.fetch():
#    fields = line.strip().split()
    if read.query_name in full_length_reads:
        fulloutfile.write(read)
	full_length += 1.0
    else:
        shortoutfile.write(read)
	short += 1.0

total = full_length + short
print sys.argv[2]
print "Total Aligned Reads:\t%d" %(int(total))
print "Total Full Length:\t%d" %(int(full_length))
print "Total Short\t%d" %(int(short))
print "%d/%d = %f %% full length" %(int(full_length),int(total),100.0*full_length/total)
print "%d/%d = %f %% short" %(int(short),int(total),100*short/total)
