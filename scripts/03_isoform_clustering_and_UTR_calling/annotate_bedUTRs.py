#!/usr/bin/env python2

## edits .bed files to reflect putative 3'UTR sites based on predicted stop codon. (Predicted by identifying largest open reading frame in clustered reads)
##
##
import sys

#"%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%i\n" %(gene, chrom, cds_strand, cds, seq, cds_start, cds_end, num_reads, num_exons, stop_codon, ",".join(stop_positions), ",".join(read_ids),retained_intron))
## Set basic parameters:

last_gene = None
last_chrom = None
last_strand = None
stop_position_dict = {}
read_to_start = {}
read_to_stop = {}
infile = open(sys.argv[1])
bedfile = open(sys.argv[2])
outfile = open(sys.argv[3],'w')
for line in infile:
    fields = line.strip().split()
    gene = fields[0]
    chrom = fields[1]
    strand = fields[2]
    cds = fields[3]
    seq = fields[4]
    cds_start = int(fields[5])
    cds_end = int(fields[6])
    num_reads = int(fields[7])
    num_exons = int(fields[8])
    start_codon = int(fields[9])
    stop_codon = int(fields[10])
    start_positions = [int(x) for x in fields[11].split(',')]
    stop_positions = [int(x) for x in fields[12].split(',')]
    
    read_ids = fields[13].split(',')
    retained_intron = fields[14] == '1'
    assert len(stop_positions) == len(read_ids)
    for read_id in read_ids:
        read_to_stop[read_id] = stop_codon
        read_to_start[read_id] = start_codon

for line in bedfile:
    fields = line.strip().split()
    newfields = fields
    read_id = fields[3]
    start = int(fields[1])
    end = int(fields[2])
    if read_id in read_to_stop:
        stop_codon = read_to_stop[read_id]
        start_codon = read_to_start[read_id]
    else:
        continue
    if fields[5] == '+':
        if start_codon > start and start_codon < end:
            newfields[6] = str(start_codon)
        else:
            newfields[6] = str(start)
        if stop_codon < end and stop_codon > start:
            newfields[7] = str(stop_codon)
        else:
            newfields[7] = str(end)
    elif fields[5] == '-':
        if stop_codon > start and stop_codon < end:
            newfields[6] = str(stop_codon)
        else:
            newfields[6] = str(start)
        if start_codon < end and start_codon > start:
            newfields[7] = str(start_codon)
        else:
            newfields[7] = str(end)
    assert int(newfields[6]) < int(newfields[7])
    
    newline = "\t".join(newfields)
    outfile.write("%s\n" %(newline))