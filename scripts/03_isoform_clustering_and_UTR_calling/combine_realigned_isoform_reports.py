#!/usr/bin/env python2
import sys

l1_in = open("../../results/isoforms/L1_isoforms.tsv")
l2_in = open("../../results/isoforms/L2_isoforms.tsv")
l3_in = open("../../results/isoforms/L3_isoforms.tsv")
l4_in = open("../../results/isoforms/L4_isoforms.tsv")
ya_in = open("../../results/isoforms/YA_isoforms.tsv")
ga_in = open("../../results/isoforms/GA_isoforms.tsv")
ml_in = open("../../results/isoforms/ML_isoforms.tsv")

all_dict = {}

for line in l1_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        introns = fields[15]
    else:
        introns = "NotSpliced"
    if gene_id in all_dict:
        if introns in all_dict[gene_id]: #Doesn't currently do 5' truncation tolerance
            num_reads = int(all_dict[gene_id][introns][7]) + int(fields[7])
            all_dict[gene_id][introns][7] = str(num_reads)
            all_dict[gene_id][introns][11] = all_dict[gene_id][introns][11] + ',' + fields[11]
            all_dict[gene_id][introns][12] = all_dict[gene_id][introns][12] + ',' + fields[12]
            all_dict[gene_id][introns][13] = all_dict[gene_id][introns][13] + ',' + fields[13]
            if len(fields[4]) > len(all_dict[gene_id][introns][4]): #replace representative read if need be
                all_dict[gene_id][introns][3] = fields[3]
                all_dict[gene_id][introns][4] = fields[4]
                all_dict[gene_id][introns][5] = fields[5]
                all_dict[gene_id][introns][6] = fields[6]
                all_dict[gene_id][introns][9] = fields[9]
                all_dict[gene_id][introns][10] = fields[10]
        else:
            all_dict[gene_id][introns] = fields
    else:
        all_dict[gene_id] = {introns:fields}
        
for line in l2_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        introns = fields[15]
    else:
        introns = "NotSpliced"
    if gene_id in all_dict:
        if introns in all_dict[gene_id]: #Doesn't currently do 5' truncation tolerance
            num_reads = int(all_dict[gene_id][introns][7]) + int(fields[7])
            all_dict[gene_id][introns][7] = str(num_reads)
            all_dict[gene_id][introns][11] = all_dict[gene_id][introns][11] + ',' + fields[11]
            all_dict[gene_id][introns][12] = all_dict[gene_id][introns][12] + ',' + fields[12]
            all_dict[gene_id][introns][13] = all_dict[gene_id][introns][13] + ',' + fields[13]
            if len(fields[4]) > len(all_dict[gene_id][introns][4]): #replace representative read if need be
                all_dict[gene_id][introns][3] = fields[3]
                all_dict[gene_id][introns][4] = fields[4]
                all_dict[gene_id][introns][5] = fields[5]
                all_dict[gene_id][introns][6] = fields[6]
                all_dict[gene_id][introns][9] = fields[9]
                all_dict[gene_id][introns][10] = fields[10]
        else:
            all_dict[gene_id][introns] = fields
    else:
        all_dict[gene_id] = {introns:fields}
        
for line in l3_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        introns = fields[15]
    else:
        introns = "NotSpliced"
    if gene_id in all_dict:
        if introns in all_dict[gene_id]: #Doesn't currently do 5' truncation tolerance
            num_reads = int(all_dict[gene_id][introns][7]) + int(fields[7])
            all_dict[gene_id][introns][7] = str(num_reads)
            all_dict[gene_id][introns][11] = all_dict[gene_id][introns][11] + ',' + fields[11]
            all_dict[gene_id][introns][12] = all_dict[gene_id][introns][12] + ',' + fields[12]
            all_dict[gene_id][introns][13] = all_dict[gene_id][introns][13] + ',' + fields[13]
            if len(fields[4]) > len(all_dict[gene_id][introns][4]): #replace representative read if need be
                all_dict[gene_id][introns][3] = fields[3]
                all_dict[gene_id][introns][4] = fields[4]
                all_dict[gene_id][introns][5] = fields[5]
                all_dict[gene_id][introns][6] = fields[6]
                all_dict[gene_id][introns][9] = fields[9]
                all_dict[gene_id][introns][10] = fields[10]
        else:
            all_dict[gene_id][introns] = fields
    else:
        all_dict[gene_id] = {introns:fields}
        
for line in l4_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        introns = fields[15]
    else:
        introns = "NotSpliced"
    if gene_id in all_dict:
        if introns in all_dict[gene_id]: #Doesn't currently do 5' truncation tolerance
            num_reads = int(all_dict[gene_id][introns][7]) + int(fields[7])
            all_dict[gene_id][introns][7] = str(num_reads)
            all_dict[gene_id][introns][11] = all_dict[gene_id][introns][11] + ',' + fields[11]
            all_dict[gene_id][introns][12] = all_dict[gene_id][introns][12] + ',' + fields[12]
            all_dict[gene_id][introns][13] = all_dict[gene_id][introns][13] + ',' + fields[13]
            if len(fields[4]) > len(all_dict[gene_id][introns][4]): #replace representative read if need be
                all_dict[gene_id][introns][3] = fields[3]
                all_dict[gene_id][introns][4] = fields[4]
                all_dict[gene_id][introns][5] = fields[5]
                all_dict[gene_id][introns][6] = fields[6]
                all_dict[gene_id][introns][9] = fields[9]
                all_dict[gene_id][introns][10] = fields[10]
        else:
            all_dict[gene_id][introns] = fields
    else:
        all_dict[gene_id] = {introns:fields}
        
for line in ya_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        introns = fields[15]
    else:
        introns = "NotSpliced"
    if gene_id in all_dict:
        if introns in all_dict[gene_id]: #Doesn't currently do 5' truncation tolerance
            num_reads = int(all_dict[gene_id][introns][7]) + int(fields[7])
            all_dict[gene_id][introns][7] = str(num_reads)
            all_dict[gene_id][introns][11] = all_dict[gene_id][introns][11] + ',' + fields[11]
            all_dict[gene_id][introns][12] = all_dict[gene_id][introns][12] + ',' + fields[12]
            all_dict[gene_id][introns][13] = all_dict[gene_id][introns][13] + ',' + fields[13]
            if len(fields[4]) > len(all_dict[gene_id][introns][4]): #replace representative read if need be
                all_dict[gene_id][introns][3] = fields[3]
                all_dict[gene_id][introns][4] = fields[4]
                all_dict[gene_id][introns][5] = fields[5]
                all_dict[gene_id][introns][6] = fields[6]
                all_dict[gene_id][introns][9] = fields[9]
                all_dict[gene_id][introns][10] = fields[10]
        else:
            all_dict[gene_id][introns] = fields
    else:
        all_dict[gene_id] = {introns:fields}
        
for line in ga_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        introns = fields[15]
    else:
        introns = "NotSpliced"
    if gene_id in all_dict:
        if introns in all_dict[gene_id]: #Doesn't currently do 5' truncation tolerance
            num_reads = int(all_dict[gene_id][introns][7]) + int(fields[7])
            all_dict[gene_id][introns][7] = str(num_reads)
            all_dict[gene_id][introns][11] = all_dict[gene_id][introns][11] + ',' + fields[11]
            all_dict[gene_id][introns][12] = all_dict[gene_id][introns][12] + ',' + fields[12]
            all_dict[gene_id][introns][13] = all_dict[gene_id][introns][13] + ',' + fields[13]
            if len(fields[4]) > len(all_dict[gene_id][introns][4]): #replace representative read if need be
                all_dict[gene_id][introns][3] = fields[3]
                all_dict[gene_id][introns][4] = fields[4]
                all_dict[gene_id][introns][5] = fields[5]
                all_dict[gene_id][introns][6] = fields[6]
                all_dict[gene_id][introns][9] = fields[9]
                all_dict[gene_id][introns][10] = fields[10]
        else:
            all_dict[gene_id][introns] = fields
    else:
        all_dict[gene_id] = {introns:fields}
        
        
for line in ml_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        introns = fields[15]
    else:
        introns = "NotSpliced"
    if gene_id in all_dict:
        if introns in all_dict[gene_id]: #Doesn't currently do 5' truncation tolerance
            num_reads = int(all_dict[gene_id][introns][7]) + int(fields[7])
            all_dict[gene_id][introns][7] = str(num_reads)
            all_dict[gene_id][introns][11] = all_dict[gene_id][introns][11] + ',' + fields[11]
            all_dict[gene_id][introns][12] = all_dict[gene_id][introns][12] + ',' + fields[12]
            all_dict[gene_id][introns][13] = all_dict[gene_id][introns][13] + ',' + fields[13]
            if len(fields[4]) > len(all_dict[gene_id][introns][4]): #replace representative read if need be
                all_dict[gene_id][introns][3] = fields[3]
                all_dict[gene_id][introns][4] = fields[4]
                all_dict[gene_id][introns][5] = fields[5]
                all_dict[gene_id][introns][6] = fields[6]
                all_dict[gene_id][introns][9] = fields[9]
                all_dict[gene_id][introns][10] = fields[10]
        else:
            all_dict[gene_id][introns] = fields
    else:
        all_dict[gene_id] = {introns:fields}
        

for gene_id in all_dict:
    for intron_chain in all_dict[gene_id]:
        print "\t".join(all_dict[gene_id][intron_chain])