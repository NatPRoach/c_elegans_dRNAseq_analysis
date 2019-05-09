#!/usr/bin/env python2
#script written by Nathan Roach nroach2@jhu.edu
#Purpose is to convert gff files in the Wormbase format to UCSC format
# ID=Gene:WBGene00000001; Name=WBGene000000001 to gene_id "WBGene000000001"; transcript_id "Y74C9A.6"
import sys

def processAttributes(raw):
    out = {}
    split1 = raw.split(';')
    for attr in split1:
        split2 = attr.split("=")
        out[split2[0]] = split2[1]
    return out

# def reformatAttributes(attr):
#     attr_list = []
#     for key in attr:
#         if key == "gene_id" or key == "transcript_id":
#             attr_list.append(key + " \"" + attr[key] +'\"')
#         else:
#             attr_list.append(key + "=" + attr[key])
#     return ";".join(attr_list)

infile = open(sys.argv[1])

cdsIDs = {}
for line in infile:
    fields = line.strip().split()
    attributes = processAttributes(fields[8])
    if fields[2] == "CDS":
        cds_id = attributes["ID"].split(":")[1]
        chrom = fields[0]
        start = int(fields[3])
        end = int(fields[4])
        strand = fields[6]
        if cds_id in cdsIDs:
            if start < cdsIDs[cds_id][1]:
                cdsIDs[cds_id][1] = start
            if end > cdsIDs[cds_id][2]:
                cdsIDs[cds_id][2] = end
        else:
            cdsIDs[cds_id] = [chrom,start,end,strand]
    else:
        continue
        
for cds in cdsIDs:
    chrom,start,end,strand = cdsIDs[cds]
    if chrom == "MtDNA":
        chrom = "chrM"
    else:
        if len(chrom) > 3:
            if chrom[0:3] != "chr":
                chrom = "chr" + chrom
        else:
            chrom = "chr" + chrom
    start = start - 1
    print "%s\t%d\t%d\t%s\t0\t%s" %(chrom,start,end,cds,strand)
    
    