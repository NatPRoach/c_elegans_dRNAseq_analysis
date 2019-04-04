#!/usr/bin/env python

import sys
import pybedtools
#import matplotlib.pyplot as pyplot
bamfilepath = sys.argv[1]
bedfilepath = sys.argv[2]
coverage_file = open(sys.argv[3])
outfile = open(sys.argv[4] + "_meta.txt",'w')

bedfile = open(bedfilepath,'r')
#Convert bed to searchable data structure:
genes = {}
duplicates = {}
last = ()
last_gene = ""
for line in bedfile:
    fields = line.split('\t')
    gene = fields[3]
    if gene not in genes:
        # exons, CD start, CD end, counts,gene_len
        genes[gene] = [[],None,None,[],0]
        duplicates[gene] = []

    if fields[7] == "exon":
        genes[gene][0].append((int(fields[1]),int(fields[2])))
        genes[gene][4] += int(fields[2]) - int(fields[1]) #update total length
        if last == (fields[1],fields[2]):
            duplicates[gene].append((int(fields[1]),int(fields[2])))
        last_gene = gene
        last = (fields[1],fields[2])
    elif fields[7] == "start_codon":
        if fields[5] == '+':
            field = 1
        else:
            field = 2
        genes[gene][1] = int(fields[field])
    elif fields[7] == "stop_codon":
        if fields[5] == '+':
            field = 1
        else:
            field = 2
        genes[gene][2] = int(fields[field])
bedfile.close()

#5'
five_prime_window = [0 for x in range(50000)]
#CDS
cds_window = [0 for x in range(100000)]
#3'
three_prime_window = [0 for x in range(50000)]

print("Done building searchable gene fields")

#a = pybedtools.BedTool(bedfilepath)
#b = pybedtools.BedTool(bamfilepath)
#c = a.coverage(b,split=True,d=True)

#for line in c:
#    print(line.strip())
    #fields = line.strip().split('\t')
    #if len(line) < 12:
    #    continue
#    if line[7] == "exon":
#        gene = line[3]
#        genes[gene][3].append(int(line[-1]))

for line in coverage_file:
    #print(line.strip())
    fields = line.strip().split('\t')
    #if len(line) < 12:
    #    continue
    if fields[7] == "exon":
        gene = fields[3]
        genes[gene][3].append(int(fields[-1]))

print("Done processing bedtools output")
gene_num = len(genes)
for gene in genes:
    geneinfo = genes[gene]

    cds_start_offset = 0
    cds_end_offset = 0
    if geneinfo[1] is not None and geneinfo[2] is not None:
        if geneinfo[1] < geneinfo[2]: #+strand
#            continue
            first_base = geneinfo[0][0][0]
            last_base = geneinfo[0][-1][1]
            for exon in geneinfo[0]:
                if geneinfo[1] > exon[1]:
                    cds_start_offset += (exon[1] - exon[0])
                elif geneinfo[1] >= exon[0]:
                    cds_start_offset += (geneinfo[1] - exon[0])
                    break
                else:
                    raise("Something went wrong")
            for exon in reversed(geneinfo[0]):
                if geneinfo[2] < exon[0]:
                    cds_end_offset += (exon[1] - exon[0])
                elif geneinfo[2] <= exon[1]:
                    cds_end_offset += (exon[1] - geneinfo[2])
                    break
                else:
                    raise("Something went wrong")
        else: # - strand
#            continue
            first_base = geneinfo[0][-1][1]
            last_base = geneinfo[0][0][0]
            geneinfo[3].reverse()
            for exon in reversed(geneinfo[0]):
                if geneinfo[1] < exon[0]:
                    cds_start_offset += (exon[1] - exon[0])
                elif geneinfo[1] <= exon[1]:
                    cds_start_offset += (exon[1] - geneinfo[1])
                    break
                else:
                    raise("Something went wrong")
            for exon in geneinfo[0]:
                if geneinfo[2] > exon[1]:
                    cds_end_offset += (exon[1] - exon[0])
                elif geneinfo[2] >= exon[0]:
                    cds_end_offset += ( geneinfo[2] - exon[0])
                    break
                else:
                    raise("Something went wrong")

        cds_end_offset = geneinfo[4] - cds_end_offset
        if cds_end_offset > cds_start_offset and cds_start_offset != 0 and cds_end_offset != geneinfo[4]:
            block = 50000.0/cds_start_offset
            for i, count in enumerate(geneinfo[3][:cds_start_offset]):
                for j in range(int(round(i*block)),int(round((i+1)*block))):
                    five_prime_window[j] += count

            block = 50000.0/(geneinfo[4]-cds_end_offset)
            for i, count in enumerate(geneinfo[3][cds_end_offset:]):
                for j in range(int(round(i*block)),int(round((i+1)*block))):
                    three_prime_window[j] += count

            block = 100000.0/(cds_end_offset - cds_start_offset)
            for i, count in enumerate(geneinfo[3][cds_start_offset:cds_end_offset]):
                for j in range(int(round(i*block)),int(round((i+1)*block))):
                    cds_window[j] += count
print("Done calculating metagene info")

data = five_prime_window + cds_window + three_prime_window
totalCounts = sum(data)
norm_data = [float(x) / float(totalCounts) for x in data]
for value in norm_data:
    outfile.write("%.64f\n" %(value))

print("Done writing data file")