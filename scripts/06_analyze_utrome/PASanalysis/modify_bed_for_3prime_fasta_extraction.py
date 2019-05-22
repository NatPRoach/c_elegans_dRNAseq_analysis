#!/usr/bin/env python2

import sys

## For now assume bed12
infile = open(sys.argv[1])
## Outputting bed6 in v2
outfile = open(sys.argv[2],'w')
window = 60 

for line in infile:
    fields = line.strip().split()
    start = int(fields[1])
    end  = int(fields[2])
    strand = fields[5]
    if end - start < window:
        difference = window - (end - start)
        if strand == '+':
            start -= difference
            end += window
        elif strand == '-':
            end += difference
            start -= window
    else:
        if strand == '+':
            end += window
        elif strand == '-':
            start -= window
    start = str(start)
    end = str(end)
    line_out = "%s\t%s\t%s\t%s\t%s\t%s\n" %(fields[0],start,end,fields[3],fields[4],fields[5])
    outfile.write(line_out)