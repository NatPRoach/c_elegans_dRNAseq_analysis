#!/usr/bin/env python2

import sys
import sets
import matplotlib
import matplotlib.pyplot as plt

def splitAttributes(attr):
    fields = attr.split(';')
    attr_dict = {}
    for field in fields:
        if len(field) == 0:
            continue
        key,val = field.strip().split(' ')
        val = val.strip('\"')
        attr_dict[key] = val
    return attr_dict

infile = open("../../references/WBcel235/Caenorhabditis_elegans.WBcel235.93.gtf")

gene_to_transcripts = {}
for line in infile:
    if line[0:2] == '#!':
        continue
    fields = line.strip().split('\t')
    chrom = fields[0]
    if chrom == "MtDNA":
        chrom = "chrM"
    else:
        chrom = "chr" + chrom
    start = int(fields[3]) - 1
    end =  int(fields[4])
    strand = fields[6]
    if fields[2] == "gene":
        attr = splitAttributes(fields[8])
        if "gene_id" in attr:
            gene_id = attr["gene_id"]
        else:
            print "something went wrong"
    if fields[2] == "transcript":
        attr = splitAttributes(fields[8])
        if "gene_id" in attr:
            gene_id = attr["gene_id"]
        else:
            print "something went wrong"
        if "transcript_id" in attr:
            tx_id = attr["transcript_id"]
        else:
            print "something went wrong"
        if "transcript_biotype" in attr:
            biotype = attr["transcript_biotype"]
        else:
            print "something went wrong"
        if biotype == "piRNA":
            continue
        if biotype == "miRNA":
            continue
        if biotype == "tRNA":
            continue
        if biotype == "snoRNA":
            continue
        if biotype == "rRNA":
            continue
        if biotype == "pseudogene":
            continue
        if biotype == "ncRNA":
            continue
        if gene_id in gene_to_transcripts:
            gene_to_transcripts[gene_id].add(tx_id)
        else:
            gene_to_transcripts[gene_id] = sets.Set([tx_id])
        

isoform_counts = []
print len(gene_to_transcripts)
for gene_id in gene_to_transcripts:
    #print gene_to_transcripts[gene_id]
    isoform_counts.append(len(gene_to_transcripts[gene_id]))

font = font = {"family":"sans-serif",
            "sans-serif":["Helvetica"],
            "size":10}
        
matplotlib.rc('font', **font)
plt.figure(num=None,figsize=(4.1667,2))
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.hist(isoform_counts,bins = [x - 0.5 for x in range(1,10)],density=True,rwidth=0.5)
plt.xticks(range(1,10))
#plt.title("annotation")
plt.xlabel("number of isoforms observed")
plt.ylabel("proportion of genes")
plt.tight_layout()
plt.savefig("../../figures/supplementals/sfigure3/sfigure3B.png")