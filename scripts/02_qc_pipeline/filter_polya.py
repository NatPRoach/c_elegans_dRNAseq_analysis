#!/usr/bin/env python2

import sys
import sets
import pysam
# insertion_limit = 20
# soft_clip_limit = 15

polya_infile = open(sys.argv[1])
bam_infile = pysam.AlignmentFile(sys.argv[2],'rb')
bam_outfile = pysam.AlignmentFile(sys.argv[3],'wb',template=bam_infile)

valid_reads = sets.Set()
for i,line in enumerate(polya_infile):
    if i == 0:
        continue
    fields = line.strip().split()
    read_id = fields[0]
    qc_tag = fields[9]
    if qc_tag == "READ_FAILED_LOAD":
        continue
    elif qc_tag == "SUFFCLIP":
        continue
    elif qc_tag == "NOREGION":
        continue
    valid_reads.add(read_id)
    
for read in bam_infile.fetch():
    if read.query_name in valid_reads:
        bam_outfile.write(read)
