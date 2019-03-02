#!/usr/bin/env python2
import sys

def testOverlap(exon,intronList):
    for intron in intronList: #intron list must be sorted
        if exon[0] <= intron[0] and exon[1] >= intron[1]:
            return True
        elif exon[1] <= intron[0]:
            return False

def testOverlap2(exon,intronList):
    vals = []
    for i,intron in enumerate(intronList): #intron list must be sorted
        if exon[0] <= intron[0] and exon[1] >= intron[1]:
            vals.append(i)
        elif exon[1] <= intron[0]:
            return vals
    return vals

def getExons(start,sizes,starts):
    exons = []
    assert len(sizes) == len(starts)
    for x in range(len(sizes)):
        exons.append((start + starts[x],start+starts[x]+sizes[x]))
    return exons

def getIntrons(strand,start,sizes,starts):
    introns = []
    assert len(sizes) == len(starts)
    if strand == '+':
        for x in range(len(sizes)-1):
            introns.append((start+starts[x]+sizes[x],start+starts[x+1]))
    elif strand == '-':
        for x in range(len(sizes)-1,0,-1):
            introns.append((start+starts[x],start+starts[x-1]+sizes[x-1]))
    return tuple(introns)

def getBlocks(sizes,starts):
    int_sizes = []
    int_starts = []
    for size in sizes.split(','):
        int_sizes.append(int(size))
    for start in starts.split(','):
        int_starts.append(int(start))
    return int_sizes, int_starts
    
    
infile = open(sys.argv[1],'r')
outfile1 = open(sys.argv[2],'w')
outfile2 = open(sys.argv[3],'w')
#### 1 - Load intron structure:
intron_annotations = {}
exon_annotations = {}
gff_in = open("/Users/nproach/Documents/NPR_Notebook/01_Scripts/c_elegans.PRJNA13758.WS265.WormBase.gff3")
for line in gff_in:
    fields = line.strip().split()
    if fields[2] == "intron":
        chrom = fields[0]
        start = int(fields[3]) - 1
        end = int(fields[4])
        strand = fields[6]
        if chrom == "MtDNA":
            chrom = "chrM"
        else:
            if len(chrom) > 3:
                if chrom[0:3] != "chr":
                    chrom = "chr" + chrom
            else:
                chrom = "chr" + chrom
        
        if (chrom,strand) in intron_annotations:
            intron_annotations[(chrom,strand)].append((start,end))
        else:
            intron_annotations[(chrom,strand)] = [(start,end)]
    if fields[2] == "exon":
        chrom = fields[0]
        start = int(fields[3]) - 1
        end = int(fields[4])
        strand = fields[6]
        if chrom == "MtDNA":
            chrom = "chrM"
        else:
            if len(chrom) > 3:
                if chrom[0:3] != "chr":
                    chrom = "chr" + chrom
            else:
                chrom = "chr" + chrom
        if (chrom,strand) in exon_annotations:
            exon_annotations[(chrom,strand)].append((start,end))
        else:
            exon_annotations[(chrom,strand)] = [(start,end)]
        

for key in intron_annotations:
    intron_annotations[key].sort()
for key in exon_annotations:
    exon_annotations[key].sort()
    if key in intron_annotations:
        for exon in exon_annotations[key]:
            vals = testOverlap2(exon,intron_annotations[key])
            vals.sort(reverse=True)
            for val in vals:
                del intron_annotations[key][val]

print "Done loading intron definitions"


for line in infile:
    fields = line.strip().split()
    chrom = fields[0]
    strand = fields[5]
    start = int(fields[1])
    block_sizes,block_starts = getBlocks(fields[10],fields[11])
    exons = getExons(start,block_sizes,block_starts)
    intron_test = False
    if (chrom,strand) in intron_annotations:
        for exon in exons:
            if testOverlap(exon,intron_annotations[(chrom,strand)]):
                intron_test = True
                break
    if intron_test:
        outfile2.write(line)
    else:
        outfile1.write(line)
    