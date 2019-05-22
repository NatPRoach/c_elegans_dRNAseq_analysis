#!/usr/bin/env python

import sys
import sets


def splitAttributes(attr):
    fields = attr.split(';')
    attr_dict = {}
    for field in fields:
        key,val = field.split('=')
        attr_dict[key] = val
    return attr_dict



in_file = open(sys.argv[1],'r')
out_file = open(sys.argv[2],'w')
reference_in = open(sys.argv[3],'r')


##Construct dictionary of donors and acceptor pairs
donors = {} # chr,strand -> set of donors on that chr strand
acceptors = {} #chr,strand -> set of acceptors on that chr strand
tx_id_to_gene = {}
donor_genes = {} #chr,strand,pos -> gene_id
acceptor_genes = {} #chr,strand,pos -> gene_id
tx_dict = {}
ambiguous_donors = sets.Set() #chr,strand,pos for donors with more than one gene assignment
ambiguous_acceptors = sets.Set() #chr,strand,pos for acceptors with more than one gene assignment
single_exon_genes = {}
tx_id_to_exons = {}
window = 100

for line in reference_in:
    fields = line.strip().split()
    chrom = fields[0]
    if chrom == "MtDNA":
        chrom = "chrM"
    else:
        chrom = "chr" + chrom
    start = int(fields[3]) - 1
    end =  int(fields[4])
    strand = fields[6]
    if fields[2] == "exon":
        tx_id = fields[8].split(':')[1]
        index = tx_id.find('.') + 1
        if tx_id in tx_id_to_exons:
            tx_id_to_exons[tx_id].append((chrom,strand,start,end))
        else:
            tx_id_to_exons[tx_id] = [(chrom,strand,start,end)]
    elif fields[2] == "gene":
        attr = splitAttributes(fields[8])
        if "Name" in attr:
            gene_id = attr["Name"]
            if "sequence_name" in attr:
                tx_id = attr["sequence_name"]
                tx_id_to_gene[tx_id] = gene_id
            else:
                print "something went wrong"
        else:
            print "something went wrong"
            


start_site_dict = {}
for tx_id in tx_id_to_exons: ##WARNING: GFF3 file must be sorted for this approach to work
    exons = tx_id_to_exons[tx_id]
    index = tx_id.find('.') + 1
    t_flag = True
    while True:
        try:
            if t_flag: #For some reason some genes have a t instead of a # immediately after the 1st '.' because this format sucks
                if tx_id[index] != 't': #Need to handle that, so this is my solution.
                    int(tx_id[index])
            else:
                int(tx_id[index])
            index += 1
            t_flag = False
        except:
            break
    tx_id2 = tx_id[:index]
    gene_id = tx_id_to_gene[tx_id2]
    if len(exons) > 1:
        chrom = exons[0][0]
        strand = exons[0][1]
        if strand == '+':
            start = exons[0][2]
        elif strand == '-':
            start = exons[-1][3]
        else:
            print "something went wrong"
    else: #Ideally we'd want a way to match single exon reads with single exon genes
        chrom,strand,start,end = exons[0]
        if strand == '-':
            tmp = start
            start = end
            end = tmp
    if (chrom,strand) not in start_site_dict:
        start_site_dict[(chrom,strand)] = sets.Set()
    for x in range(start - window, start + window + 1):
        start_site_dict[(chrom,strand)].add(x)

for line in in_file:
    fields = line.strip().split()
    chrom = fields[0]
    strand = fields[5]
    if strand == '+':
        start = int(fields[1])
    elif strand == '-':
        start = int(fields[2])
    else:
        print "something went wrong"
    if (chrom,strand) in start_site_dict:
        if start in start_site_dict[(chrom,strand)]:
            out_file.write(line)

