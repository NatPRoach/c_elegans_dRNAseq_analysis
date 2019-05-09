#!/usr/bin/env python2

import sys
import pysam
import sets
import regex as re
import matplotlib.pyplot as plt
from Bio import Seq

def splitORF(input_seq):
    startP = re.compile('ATG')
    nuc = input_seq.replace('\n','')
    longest = (0,)
    for m in startP.finditer(nuc, overlapped=True):
        mod = len(Seq.Seq(nuc)[m.start():]) % 3
        if mod != 0:
            pro = Seq.Seq(nuc)[m.start():-mod].translate(to_stop=True)
        else:
            pro = Seq.Seq(nuc)[m.start():].translate(to_stop=True)
        if len(pro) > longest[0]:
            flag = True
            stop = nuc[m.start()+len(pro)*3: m.start()+len(pro)*3 +3 ]
            if len(stop) != 3:
                flag = False
            if stop != "TAG" and stop != "TAA" and stop != "TGA":
                flag = False
            if flag:
                longest = (len(pro), 
                           m.start(), 
                           m.start()+len(pro)*3+3,
                           str(pro),
                           len(nuc[m.start()+len(pro)*3+3:]))# everything after stop codon
    return longest

def equalORFs(orf1,orf2):
    if len(orf2) > len(orf1): #len(orf1) >= len(orf2)
        tmp = orf1
        orf1 = orf2
        orf2 = tmp
    return orf1[-len(orf2):] == orf2

## Removing 5' tolerance
# def equalIntrons(introns1,introns2):
#     if len(introns2) > len(introns1):
#         tmp = introns1
#         introns1 = introns2
#         introns2 = tmp
#     return introns1[-len(introns2):] == introns2
def equalIntrons(introns1,introns2):
    return introns1 == introns2



def getStopCodon(strand,start, block_sizes,block_starts,utr_len):
    if strand == '+':
        for x in range(len(block_starts)):
            if utr_len > block_sizes[len(block_starts)-1-x]:
                utr_len -= block_sizes[len(block_starts)-1-x]
            else:
                return start + block_starts[len(block_starts)-1-x] + (block_sizes[len(block_starts)-1-x] - utr_len)
    elif strand == '-':
        for x in range(len(block_starts)):
            if utr_len > block_sizes[x]:
                utr_len -= block_sizes[x]
            else:
                return start + block_starts[x] + utr_len
def getStartCodon(strand,start,block_sizes,block_starts,cds_start):
    if strand == '+':
        for x in range(len(block_starts)):
            if cds_start > block_sizes[x]:
                cds_start -= block_sizes[x]
            else:
                return start + block_starts[x] + cds_start
    elif strand == '-':
        for x in range(len(block_starts)):
            if cds_start > block_sizes[len(block_starts)-1-x]:
                cds_start -= block_sizes[len(block_starts)-1-x]
            else:
                return start + block_starts[len(block_starts)-1-x] + (block_sizes[len(block_starts)-1-x] - cds_start)
                
def getBlocks(sizes,starts):
    int_sizes = []
    int_starts = []
    for size in sizes.split(','):
        int_sizes.append(int(size))
    for start in starts.split(','):
        int_starts.append(int(start))
    return int_sizes, int_starts

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

# def testOverlap(exon,intronList): #Binary search for intron overlap
#       #Somthing wrong with the logic to this, but I'm not sure what. It doesn't yield the same overlaps as testOverlap2
#     while True:
#         if len(intronList) == 0:
#             return False
#         index = len(intronList) // 2
#         intron = intronList[index]
#         if exon[0] <= intron[0] and exon[1] >= intron[1]:
#             return True
#         if index == 0:
#             return False
#         if exon < intron:
#             intronList = intronList[:index]
#         elif exon > intron:
#             intronList = intronList[index+1:]

window = 50
cds_threshold = 1

#### 1 - Load tx_id to gene_id
fa_in = open("../../references/WS265/c_elegans.PRJNA13758.WS265.mRNA_transcripts.fa")
tx_to_gene_id = {}
for line in fa_in:
    if line[0] == ">":
        fields = line.strip("\n>").split()
        tx_id = fields[0]
        gene = fields[1].split("=")[1]
        tx_to_gene_id[tx_id] = gene
print "Done loading tx_id to gene"
        
#### 2 - for each read, get the assigned gene
read_to_gene_id = {}
counter = 0
for line in open(sys.argv[1],'r'):
    fields = line.strip().split()
    read_id = fields[0]
    gene_id = fields[1]
    read_to_gene_id[read_id] = gene_id
# if sys.argv[1][-3:] == "sam":
#     sam_in = pysam.AlignmentFile(sys.argv[1],'r')
# elif sys.argv[1][-3:] == "bam":
#     sam_in = pysam.AlignmentFile(sys.argv[1],'rb')
# for read in sam_in.fetch():
#     if read.reference_id != -1:
#         if read.reference_name in tx_to_gene_id:
#             read_to_gene_id[read.query_name] = tx_to_gene_id[read.reference_name]
print "Done loading read_id to gene"

#### 3 - Load intron structure:
# intron_annotations = {}
# exon_annotations = {}
# gff_in = open("/Users/nproach/Documents/NPR_Notebook/01_Scripts/c_elegans.PRJNA13758.WS265.WormBase.gff3")
# for line in gff_in:
#     fields = line.strip().split()
#     if fields[2] == "intron":
#         chrom = fields[0]
#         start = int(fields[3]) - 1
#         end = int(fields[4])
#         strand = fields[6]
#         if chrom == "MtDNA":
#             chrom = "chrM"
#         else:
#             if len(chrom) > 3:
#                 if chrom[0:3] != "chr":
#                     chrom = "chr" + chrom
#             else:
#                 chrom = "chr" + chrom
#
#         if (chrom,strand) in intron_annotations:
#             intron_annotations[(chrom,strand)].append((start,end))
#         else:
#             intron_annotations[(chrom,strand)] = [(start,end)]
#     if fields[2] == "exon":
#         chrom = fields[0]
#         start = int(fields[3]) - 1
#         end = int(fields[4])
#         strand = fields[6]
#         if chrom == "MtDNA":
#             chrom = "chrM"
#         else:
#             if len(chrom) > 3:
#                 if chrom[0:3] != "chr":
#                     chrom = "chr" + chrom
#             else:
#                 chrom = "chr" + chrom
#         if (chrom,strand) in exon_annotations:
#             exon_annotations[(chrom,strand)].append((start,end))
#         else:
#             exon_annotations[(chrom,strand)] = [(start,end)]
#
#
# for key in intron_annotations:
#     intron_annotations[key].sort()
# for key in exon_annotations:
#     exon_annotations[key].sort()
#     if key in intron_annotations:
#         for exon in exon_annotations[key]:
#             vals = testOverlap2(exon,intron_annotations[key])
#             vals.sort(reverse=True)
#             for val in vals:
#                 del intron_annotations[key][val]
#
# print "Done loading intron definitions"

in_file  = open(sys.argv[2],'r')
out_file = open(sys.argv[3],'w')
# out_file.write("gene_id\tfull_sequence\tCDS\tCDS_start\tCDS_end\t#ofexons\n")
gene_to_isoforms = {}
gene_to_stop_codons = {}
counter = 0
### 4 - create a structure that groups read info by their associated gene
gene_to_reads = {}
for line in in_file:
    fields = line.strip().split()
    if fields[3] in read_to_gene_id:
        gene_id = read_to_gene_id[fields[3]]
    else:
        gene_id = '-'
    # if int(fields[4]) < 60:
    #     continue
    if gene_id in gene_to_reads:
        gene_to_reads[gene_id].append(fields)
    else:
        gene_to_reads[gene_id] = [fields]

exon_counts = []
gene_counter = 0
isoform_counter = 0
intron_retention_counter = 0 
total_count = 0
print "Done loading gene to reads"

### 5 - for each gene, define isoforms
for gene in gene_to_reads:
    gene_increment_flag = False
    introns_to_read_info = {}
    for read in gene_to_reads[gene]:
        ## fetch the intron structure of a read
        start = int(read[1])
        strand = read[5]
        block_sizes,block_starts = getBlocks(read[10],read[11])
        introns = getIntrons(strand,start,block_sizes,block_starts)
        ## if that intron structure has been seen already for this gene, add read info to list, else make new list
        if introns in introns_to_read_info:
            introns_to_read_info[introns].append(read)
        else:
            introns_to_read_info[introns] = [read]
    ## create a list of intron structures, and sort based on back to front ordering of items.
    introns_list = []
    for introns in introns_to_read_info:
        introns_list.append(introns)
    introns_list.sort(key = lambda x: x[::-1])
    
    cds_counter = 1
    ## combine intron structures that differ only by n terminal truncation
    ## NOTE: Creates problem when structure is intron retained at 5' most splice site.
    ## need to change the logic?
    new_introns_to_read_info = {}
    temp = [0]
    for x in range(len(introns_list) - 1):
        if not equalIntrons(introns_list[x+1],introns_list[x]):
            cds_counter += 1 # keep track of how many unique isoforms we're identifying
            new_introns = introns_list[temp[-1]]
            for val in temp:
                introns = introns_list[val]
                if new_introns in  new_introns_to_read_info:
                    for read in introns_to_read_info[introns]:
                        new_introns_to_read_info[new_introns].append(read)
                else:
                    new_introns_to_read_info[new_introns] = introns_to_read_info[introns]
            temp = [x+1]
        else:
            temp.append(x+1)
    new_introns = introns_list[temp[-1]]
    for val in temp:
        introns = introns_list[val]
        if new_introns in  new_introns_to_read_info:
            for read in introns_to_read_info[introns]:
                new_introns_to_read_info[new_introns].append(read)
        else:
            new_introns_to_read_info[new_introns] = introns_to_read_info[introns]
    
    ### For each intron structure, fetch statistics that summarize the isoform structure
    ### E.g. exon number, stop positions of supporting reads, read ids of supporting reads
    ### etc. 
    total_count += cds_counter # keep track of how many isoforms ID'd across all genes
    for introns in new_introns_to_read_info:
        seq = ""
        cds = None
        cds_start = None
        cds_end = None
        exons = None
        num_exons = None
        stop_codon = None
        chrom = None
        cds_strand = None 
        retained_intron = None
        num_reads = len(new_introns_to_read_info[introns])
        start_positions = []
        stop_positions = []
        read_ids = []
        print_flag = True
        for read in new_introns_to_read_info[introns]:
            chrom = read[0]
            start = int(read[1])
            end = int(read[2])
            read_id = read[3]
            strand = read[5]
            block_sizes,block_starts = getBlocks(read[10],read[11])
            nuc = read[12].upper()
            split = splitORF(nuc)
            if len(split) == 1:
                continue 
            utr_len = split[4]
            
            # Grab info from the longest of the reads NOTE: Again creates a problem when 5' most splice site is retained
            # As that will often be the longest read, and the isoform will be reported as intron retained even if only subset of reads
            # are intron retained. Need a way of ID-ing that.
            
            temp_exons = getExons(start,block_sizes,block_starts)
            intron_test = False
            ## Commenting this out for now, since we're filtering out intron retention at an earlier step. 
            # for exon in temp_exons:
            #     if (chrom,strand) in intron_annotations:
            #         if testOverlap(exon,intron_annotations[(chrom,strand)]):
            #             intron_test = True
            #             break


            # if gene == "WBGene00006439" or gene == "WBGene00021043" or gene == "WBGene00008504":
#                 print gene
#                 print temp_exons
#                 #print intron_annotations[chrom,strand]
#                 print intron_test
            if retained_intron is None: # If uninitialized, assign no matter what
                seq = nuc
                cds_start = split[1]
                cds_end = split[2]
                cds = split[3]
                num_exons = int(read[9])
                exons = temp_exons
                retained_intron = intron_test
                stop_codon = getStopCodon(strand,start,block_sizes,block_starts,utr_len)
                start_codon = getStartCodon(strand,start,block_sizes,block_starts,cds_start)
            elif not intron_test and retained_intron: # If initialized as retained, and non retained candidate; assign no matter what
                seq = nuc
                cds_start = split[1]
                cds_end = split[2]
                cds = split[3]
                num_exons = int(read[9])
                exons = temp_exons
                retained_intron = intron_test
                stop_codon = getStopCodon(strand,start,block_sizes,block_starts,utr_len)
                start_codon = getStartCodon(strand,start,block_sizes,block_starts,cds_start)
            elif retained_intron and intron_test: # If initialized as retained and retained candidate; check length
                if len(nuc) > len(seq):
                    seq = nuc
                    cds_start = split[1]
                    cds_end = split[2]
                    cds = split[3]
                    num_exons = int(read[9])
                    exons = temp_exons
                    stop_codon = getStopCodon(strand,start,block_sizes,block_starts,utr_len)
                    start_codon = getStartCodon(strand,start,block_sizes,block_starts,cds_start)
            elif not retained_intron and not intron_test: #If initialized as not retained and not retained, check length
                if len(nuc) > len(seq):
                    seq = nuc
                    cds_start = split[1]
                    cds_end = split[2]
                    cds = split[3]
                    num_exons = int(read[9])
                    exons = temp_exons
                    stop_codon = getStopCodon(strand,start,block_sizes,block_starts,utr_len)
                    start_codon = getStartCodon(strand,start,block_sizes,block_starts,cds_start)
            # else: if initialized as not retained and retained candidate, skip no matter what.
                
            if chrom is not None: #Throw out everything with ambiguous mapping
                if chrom != read[0]:
                    print_flag = False
                    continue
            else:
                chrom = read[0]
            if cds_strand is not None: #Throw out everything with ambiguous mapping
                if cds_strand != strand:
                    print_flag = False
                    continue
            else:
                cds_strand = strand
            # if stop_codon is not None: #Throw out everything with ambiguous mapping
            #     if stop_codon != getStopCodon(strand,start,block_sizes,block_starts,utr_len):
            #         print_flag = False
            #         continue
            # else:
            #     stop_codon = getStopCodon(strand,start,block_sizes,block_starts,utr_len)
            #
            read_ids.append(read_id)
            
            if strand == "+":
                stop_positions.append(str(end))
                start_positions.append(str(start))
            elif strand == "-":
                stop_positions.append(str(start))
                start_positions.append(str(end))
        if print_flag and len(read_ids) != 0:
            isoform_counter += 1
            intron_list = []
            for intron in introns:
                intron_list.append(str(intron[0]) + ',' + str(intron[1]))
            intron_string = ";".join(intron_list)
            out_file.write("%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%i\t%s\n" %(gene, chrom, cds_strand, cds, seq, cds_start, cds_end, num_reads, num_exons, start_codon, stop_codon, ",".join(start_positions), ",".join(stop_positions), ",".join(read_ids),retained_intron,intron_string))
            if retained_intron:
                intron_retention_counter  += 1
            exon_counts.append(num_exons)
            gene_increment_flag = True
    if gene_increment_flag:
        gene_counter += 1


print "%d CDS isoforms over %d genes" %(isoform_counter, gene_counter)
# print "%d CDS isoforms over %d genes" %(total_count, gene_counter)
print "%d CDS isoforms had retained introns" %(intron_retention_counter)
#print "%d UTR isoforms over %d genes" %(utr_count, gene_count)