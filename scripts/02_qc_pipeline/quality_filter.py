#!/usr/bin/env python2

import sys
import pysam

infile = pysam.AlignmentFile(sys.argv[1],'rb')
outfile = pysam.AlignmentFile(sys.argv[2],'wb',template=infile)

insertion_limit = 20
soft_clip_limit = 20

for read in infile.fetch():
    cigar_tuples = read.cigartuples
    if cigar_tuples is not None:
        print_flag = True
        if read.is_reverse:
            for op, length in cigar_tuples:
                if op == 5:
                    continue
                elif op == 4:
                    if length >= soft_clip_limit:
                        print_flag = False
                else:
                    break
        else:
            for op, length in reversed(cigar_tuples):
                if op == 5:
                    continue
                elif op == 4:
                    if length >= soft_clip_limit:
                        print_flag = False
                else:
                    break
            
        for op,length in cigar_tuples:
            if op == 1 and length >= insertion_limit:
                print_flag = False
        
        if print_flag:
            outfile.write(read)

