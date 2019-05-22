#!/usr/bin/env python2

import sys
import sets
import rpy2.robjects.numpy2ri
from rpy2.robjects.packages import importr
import numpy as np
import scipy.stats as sc

def load_read_id_to_length():
    read_id_to_length = {}### For each read in the polyA file assign its polyA tail length in a dictionary
    tmppolyafile = open("../../data/L1/bio1/tech1/analysis/ce11_gen_L1_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L1/bio1/tech2/analysis/ce11_gen_L1_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L2/bio1/tech1/analysis/ce11_gen_L2_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L2/bio1/tech2/analysis/ce11_gen_L2_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])            

    tmppolyafile = open("../../data/L3/bio1/tech1/analysis/ce11_gen_L3_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L3/bio1/tech2/analysis/ce11_gen_L3_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L4/bio1/tech1/analysis/ce11_gen_L4_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L4/bio1/tech2/analysis/ce11_gen_L4_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])
            
    tmppolyafile = open("../../data/young_adult/bio1/tech1/analysis/ce11_gen_young_adult_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/young_adult/bio1/tech2/analysis/ce11_gen_young_adult_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/adult/bio1/tech1/analysis/ce11_gen_adult_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/adult/bio1/tech2/analysis/ce11_gen_adult_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])
            
    tmppolyafile = open("../../data/male/bio1/tech1/analysis/ce11_gen_male_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/male/bio1/tech2/analysis/ce11_gen_male_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])
    return read_id_to_length

def getSignificances(utr_id_file,outfile,outfile2):
    gene_to_ids = {}
    gene_to_utrs = {}
    id_to_utr = {}
    utr_to_lengths = {}
    read_id_to_length = load_read_id_to_length()
    outfile2.write("gene\tutr\tlength\n")
    for line in utr_id_file:
        fields = line.strip().split()
        read_id = fields[0]
        if read_id not in read_id_to_length:
            continue
        utr_id = fields[1]
        gene_id = fields[1].split("-cluster")[0]
        if gene_id in gene_to_ids:
            gene_to_ids[gene_id].append(read_id)
        else:
            gene_to_ids[gene_id] = [read_id]
        if gene_id not in gene_to_utrs:
            gene_to_utrs[gene_id] = sets.Set()
    
        gene_to_utrs[gene_id].add(utr_id)
        if utr_id in utr_to_lengths:
            utr_to_lengths[utr_id].append(read_id_to_length[read_id])
        else:
            utr_to_lengths[utr_id] = [read_id_to_length[read_id]]
        id_to_utr[read_id] = utr_id

    for gene_id in gene_to_ids:
        uniq_utrs = gene_to_utrs[gene_id]
        utr_list = []
        utr_idx = {}
        for i, utr in enumerate(uniq_utrs):
            utr_list.append(utr)
            utr_idx[utr] = i
        utr_num = len(utr_list)
        for utr in utr_list:
            for length in utr_to_lengths[utr]:
                outfile2.write("%s\t%s\t%f\n" %(gene_id,utr,length))
        if utr_num == 2:
            k, pval = sc.kruskal(utr_to_lengths[utr_list[0]],utr_to_lengths[utr_list[1]])
        elif utr_num == 3:
            k, pval = sc.kruskal(utr_to_lengths[utr_list[0]],utr_to_lengths[utr_list[1]],utr_to_lengths[utr_list[2]])
        elif utr_num == 4:
            k, pval = sc.kruskal(utr_to_lengths[utr_list[0]],utr_to_lengths[utr_list[1]],utr_to_lengths[utr_list[2]],utr_to_lengths[utr_list[3]])
        elif utr_num == 5:
            k, pval = sc.kruskal(utr_to_lengths[utr_list[0]], utr_to_lengths[utr_list[1]], utr_to_lengths[utr_list[2]], utr_to_lengths[utr_list[3]], utr_to_lengths[utr_list[4]])
        elif utr_num == 6:
            k, pval = sc.kruskal(utr_to_lengths[utr_list[0]], utr_to_lengths[utr_list[1]], utr_to_lengths[utr_list[2]], utr_to_lengths[utr_list[3]], utr_to_lengths[utr_list[4]], utr_to_lengths[utr_list[5]])
        elif utr_num == 7:
            k, pval = sc.kruskal(utr_to_lengths[utr_list[0]], utr_to_lengths[utr_list[1]], utr_to_lengths[utr_list[2]], utr_to_lengths[utr_list[3]], utr_to_lengths[utr_list[4]], utr_to_lengths[utr_list[5]], utr_to_lengths[utr_list[6]])
        elif utr_num == 8:
            k, pval = sc.kruskal(utr_to_lengths[utr_list[0]], utr_to_lengths[utr_list[1]], utr_to_lengths[utr_list[2]], utr_to_lengths[utr_list[3]], utr_to_lengths[utr_list[4]], utr_to_lengths[utr_list[5]], utr_to_lengths[utr_list[6]], utr_to_lengths[utr_list[7]])
        else:
            continue

        outfile.write("%s\t%s\t%e\n" %(gene_id,",".join(utr_list),pval))

# getSignificances(open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/assignments/L1_utrs.tsv"),open("l1.utr_polya_sig.txt",'w'),open("l1.utr_polya_lengths.txt",'w'))
# getSignificances(open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/assignments/L2_utrs.tsv"),open("l2.utr_polya_sig.txt",'w'),open("l2.utr_polya_lengths.txt",'w'))
# getSignificances(open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/assignments/L3_utrs.tsv"),open("l3.utr_polya_sig.txt",'w'),open("l3.utr_polya_lengths.txt",'w'))
# getSignificances(open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/assignments/L4_utrs.tsv"),open("l4.utr_polya_sig.txt",'w'),open("l4.utr_polya_lengths.txt",'w'))
# getSignificances(open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/assignments/YA_utrs.tsv"),open("ya.utr_polya_sig.txt",'w'),open("ya.utr_polya_lengths.txt",'w'))
# getSignificances(open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/assignments/GA_utrs.tsv"),open("ga.utr_polya_sig.txt",'w'),open("ga.utr_polya_lengths.txt",'w'))
# getSignificances(open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/assignments/ML_utrs.tsv"),open("ml.utr_polya_sig.txt",'w'),open("ml.utr_polya_lengths.txt",'w'))
getSignificances(open("../../results/utrs/assignments/all_isoforms_utrs.tsv"),open("../../results/scratch/polya/all.utr_polya_sig.txt",'w'),open("../../results/scratch/polya/all.utr_polya_lengths.txt",'w'))