#!/usr/bin/env python2

import sets

assignments_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/analyzeUTRome/PASanalysis/assignmentFiles/all_isoforms.PASassignments.txt")
my_utrs = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/all_isoforms_utrs.bed",'r')
mangone_utrs = open("/Users/nproach/Documents/NPR_Notebook/00_Data/references/utrs/mangone/ce11_converted.stop_codon_converted.bed",'r')
bartel_utrs = open("/Users/nproach/Documents/NPR_Notebook/00_Data/references/utrs/bartel/ce11_converted.shifted.bed",'r')
# wormbase_utrs = open("/Users/nproach/Documents/NPR_Notebook/00_Data/references/utrs/wormbase/wormbase_utrs.bed",'r')
outfile = open("mangone_overlapping.bed",'w')
#outfile = open("utr.overlap.matrix",'w')
#outfile2 = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/novel/all_novel_utrs.bed",'w')
window = 10

# cluster_to_assignment = {}
# for line in assignments_in:
#     fields = line.strip().split()
#     cluster_id = fields[0]
#     pas_type = fields[1]
#     if pas_type != "AATAAA" and pas_type != "noPAS":
#         pas_type = "altPAS"
#     cluster_to_assignment[cluster_id] = pas_type
#

my_utr_dict = {}
recall_dict = {}
my_bed_dict = {}
for line in my_utrs:
    fields = line.strip().split()
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    cluster_id = fields[3]
    strand = fields[5]
    if strand == '+':
        id_tuple = (chrom,strand,start)
        id_tuple2 = (chrom,strand,start,end)
        cleavage_site = end
    elif strand == '-':
        id_tuple = (chrom,strand,end)
        id_tuple2 = (chrom,strand,end,start)
        cleavage_site = start
    else:
        print "Something went wrong"
    if id_tuple in my_utr_dict:
        my_utr_dict[id_tuple].append(cleavage_site)
    else:
        my_utr_dict[id_tuple] = [cleavage_site]
    if id_tuple2 in recall_dict:
        recall_dict[id_tuple2].append(cluster_id)
    else:
        recall_dict[id_tuple2] = [cluster_id]
    my_bed_dict[cluster_id] = line
    
mangone_utr_dict = {}
for line in mangone_utrs:
    fields = line.strip().split()
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    strand = fields[5]
    if strand == '+':
        id_tuple = (chrom,strand,start)
        cleavage_site = end
    elif strand == '-':
        id_tuple = (chrom,strand,end)
        cleavage_site = start
    else:
        print "Something went wrong"
    if id_tuple in mangone_utr_dict:
        mangone_utr_dict[id_tuple].append(cleavage_site)
    else:
        mangone_utr_dict[id_tuple] = [cleavage_site]
    
bartel_utr_dict = {}
for line in bartel_utrs:
    fields = line.strip().split()
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    strand = fields[5]
    if strand == '+':
        id_tuple = (chrom,strand,start)
        cleavage_site = end
    elif strand == '-':
        id_tuple = (chrom,strand,end)
        cleavage_site = start
    else:
        print "Something went wrong"
    if id_tuple in bartel_utr_dict:
        bartel_utr_dict[id_tuple].append(cleavage_site)
    else:
        bartel_utr_dict[id_tuple] = [cleavage_site]


total_set2 = sets.Set()
for chrom,strand,stop_codon in mangone_utr_dict:
    for cleavage_site in mangone_utr_dict[(chrom,strand,stop_codon)]:
        total_set2.add((chrom,strand,stop_codon,cleavage_site))
# for chrom,strand,stop_codon in bartel_utr_dict:
#     for cleavage_site in bartel_utr_dict[(chrom,strand,stop_codon)]:
#         total_set2.add((chrom,strand,stop_codon,cleavage_site))
# for chrom,strand,stop_codon in wormbase_utr_dict:
#     for cleavage_site in wormbase_utr_dict[(chrom,strand,stop_codon)]:
#         total_set2.add((chrom,strand,stop_codon,cleavage_site))

for chrom,strand,stop_codon in my_utr_dict:
    for cleavage_site in my_utr_dict[(chrom,strand,stop_codon)]:
        if (chrom,strand,stop_codon,cleavage_site) in total_set2:
            continue
            #wormbase_set.add((chrom,strand,stop_codon,cleavage_site))
        else:
            flag = True
            for x in range(cleavage_site - window, cleavage_site + window + 1): # is it within a window?
                # if stop_codon == 11420742:
                #     print "here1" , x
                #     if x == 11420485:
                #         print "here"
                if (chrom,strand,stop_codon,x) in total_set2:
                    t = (chrom,strand,stop_codon,cleavage_site)
                    for cluster_id in recall_dict[t]:
                        outfile.write(my_bed_dict[cluster_id])
                    flag = False
                    break
            # if flag:
            #     t = (chrom,strand,stop_codon,cleavage_site)
            #     for cluster_id in recall_dict[t]:
            #         outfile2.write(my_bed_dict[cluster_id])