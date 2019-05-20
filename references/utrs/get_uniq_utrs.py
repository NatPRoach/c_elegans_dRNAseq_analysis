#!/usr/bin/env python
import sets

infile = open("mangone_utrs.bed")
outfile = open("mangone_uniq_utrs.bed",'w')

uniq = sets.Set()
for line in infile:
    fields = line.strip().split()
    chrom = fields[0]
    start = fields[1]
    end = fields[2]
    strand = fields[5]
    if (chrom,start,end,strand) in uniq:
        continue
    else:
        uniq.add((chrom,start,end,strand))
        outfile.write(line)

infile.close()
outfile.close()

infile = open("jan_utrs.bed")
outfile = open("jan_uniq_utrs.bed",'w')

uniq = sets.Set()
for line in infile:
    fields = line.strip().split()
    chrom = fields[0]
    start = fields[1]
    end = fields[2]
    strand = fields[5]
    if (chrom,start,end,strand) in uniq:
        continue
    else:
        uniq.add((chrom,start,end,strand))
        outfile.write(line)
