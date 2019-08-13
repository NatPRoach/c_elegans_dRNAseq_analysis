#!/usr/bin/env python2

import sys
import sets

def splitAttributes(attr):
    fields = attr.split(';')
    attr_dict = {}
    for field in fields:
        key,val = field.split('=')
        attr_dict[key] = val
    return attr_dict

def getIntrons(fields):
    strand = fields[5]
    
    introns = []
    
    start = int(fields[1])
    blockSizes  = [int(x) for x in fields[10].split(',')]
    blockStarts = [int(x) for x in fields[11].split(',')]
    
    if len(blockSizes) > 1:
        for x in range(len(blockSizes)-1):
            size1 = blockSizes[x]
            size2 = blockSizes[x+1]
            start1 = start + blockStarts[x]
            end1 = start1 + size1
            start2 = start + blockStarts[x+1]
            end2 = start2 + size2
            if strand == '+':
                introns.append((end1,start2))
            elif strand == '-':
                introns.append((start2,end1))
    if strand == '-':
        introns.reverse()
    return tuple(introns)

### 1 - get annotation info, the intron chain for each transcript isoform
ingff = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/c_elegans.PRJNA13758.WS265.WormBase.exons.chr.gff3")
tx_id_to_exons = {}
for line in ingff:
    fields = line.strip().split()
    chrom = fields[0]
    if chrom == "MtDNA":
        chrom = "chrM"
    else:
        if len(chrom) < 3:
            chrom = "chr" + chrom
        elif chrom[:3] != "chr":
            chrom = "chr" + chrom
    start = int(fields[3]) - 1
    end =  int(fields[4])
    strand = fields[6]
    if fields[2] == "exon":
        tx_id = fields[8].split(':')[1]
        if tx_id in tx_id_to_exons:
            tx_id_to_exons[tx_id].append((chrom,strand,start,end))
        else:
            tx_id_to_exons[tx_id] = [(chrom,strand,start,end)]
    

tx_id_to_introns = {}
introns_to_tx_id = {}
intron_set = sets.Set()
for tx_id in tx_id_to_exons: ##WARNING: GFF3 file must be sorted for this approach to work
    exons = tx_id_to_exons[tx_id]
    if len(exons) > 1:
        introns = []
        chrom = None
        strand = None
        for x in range(len(exons)-1):
            exon1 = exons[x]
            exon2 = exons[x+1]
            assert exon1[0] == exon2[0]
            chrom = exon1[0]
            assert exon1[1] == exon2[1]
            strand = exon1[1]
            if strand == '+':
                donor = exon1[3]
                acceptor = exon2[2]
            elif strand == '-':
                donor = exon2[2]
                acceptor = exon1[3]
            introns.append((donor,acceptor))
        if strand == '+':
            introns = tuple(introns)
        elif strand == '-':
            introns.reverse()
            introns = tuple(introns)
        # #store intron structure, and all 5' truncations of intron structure and link to the corresponding tx_id
        # for x in range(len(introns)):
        #     intron_subset = introns[x:]
        #     if intron_subset in introns_to_tx_id:
        #         introns_to_tx_id[intron_subset].append(tx_id)
        #     else:
        #         introns_to_tx_id[intron_subset] = [tx_id]
        #     intron_set.add(introns[x])
        for x in range(1,len(introns)):
            intron_subset = introns[:len(introns) - x]
            intron_set.add(intron_subset)

infile = open(sys.argv[1])
total = 0
count = 0
previous_existing = 0
unambiguous = 0
new_introns = 0
for line in infile:
    fields = line.strip().split()
    introns = getIntrons(fields)
    if introns in intron_set:
        count += 1
        print line.strip()
    total += 1
    # if len(fields) == 16:
    #     if fields[14] != '1': # Ignore intron retained reads
    #         intron_string = fields[15]
    #         introns = []
    #         for pair in intron_string.split(';'):
    #             split_pair = pair.split(',')
    #             if len(split_pair) == 2:
    #                 intron = (int(split_pair[0]),int(split_pair[1]))
    #                 if intron not in intron_set:
    #                     new_introns += 1
    #                 introns.append(intron)
    #         introns = tuple(introns)
    #         if introns in introns_to_tx_id:
    #             previous_existing += 1
    #             if len(introns_to_tx_id[introns]) == 1:
    #                 unambiguous +=1
    #         else:
    #             print line.strip() # Report any novel isoforms by printing them to stdout
    #         total += 1
    # else:
    #     unambiguous += 1
    #     previous_existing += 1
    #     total += 1
# print count
# print total
sys.stderr.write("new introns: %d\n" %(new_introns))