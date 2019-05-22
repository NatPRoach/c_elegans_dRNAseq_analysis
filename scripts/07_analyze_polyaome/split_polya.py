#!/usr/bin/env python2

import sys
import numpy
import scipy.stats

txoutfile2 =   open("../../results/scratch/polya/fulls/all.full.polya"  ,'w')

## Initialize dictionaries used to store data
tx_id_to_lengths = {} #for each tx_id (# of reads, [polyA lengths])
gene_id_to_lengths = {} # for each gene_id (# of reads, [polyA lengths])
tx_id_to_gene_id = {} #for each Wormbase Tx ID -> Wormbase gene_id
gene_id_to_tx_id = {} #for each Wormbase geneID -> [Wormbase Tx IDs]
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

mRNA_tx_fasta = open("../../references/WS265/c_elegans.PRJNA13758.WS265.mRNA_transcripts.fa")
for line in mRNA_tx_fasta:
    if line[0] != '>':
        continue
    fields = line.strip().split()
    tx_id = fields[0].strip(">")
    gene_id = fields[1].split("=")[1]
    tx_id_to_gene_id[tx_id] = gene_id
    if gene_id in gene_id_to_tx_id:
        gene_id_to_tx_id[gene_id].append(tx_id)
    else:
        gene_id_to_tx_id[gene_id] = [tx_id]


tmp_tx_id_file = open("../../results/correctionLogs/l1.gene.txt")
for line in tmp_tx_id_file:
    fields = line.strip().split()
    if len(fields) != 4:
        continue
    read_id = fields[0]
    gene_id = fields[1]
    tx_ids = fields[3].split(',')
    if len(tx_ids) > 1:
        continue
    tx_id = tx_ids[0]
    if read_id in read_id_to_length:
        if tx_id in tx_id_to_lengths: #No need to initialize
            tx_id_to_lengths[tx_id][0] += 1
            tx_id_to_lengths[tx_id][1].append(read_id_to_length[read_id])
            gene_id_to_lengths[gene_id][0] += 1
            gene_id_to_lengths[gene_id].append(read_id_to_length[read_id])
        elif gene_id in gene_id_to_lengths: #other isoform seen?
            tx_id_to_lengths[tx_id] = [1, [read_id_to_length[read_id]]]
            gene_id_to_lengths[gene_id][0] += 1
            gene_id_to_lengths[gene_id].append(read_id_to_length[read_id])
        else: #initialize both
            tx_id_to_lengths[tx_id] = [1, [read_id_to_length[read_id]]]
            gene_id_to_lengths[gene_id] = [1, [read_id_to_length[read_id]]]


tmp_tx_id_file = open("../../results/correctionLogs/l2.gene.txt")
for line in tmp_tx_id_file:
    fields = line.strip().split()
    if len(fields) != 4:
        continue
    read_id = fields[0]
    gene_id = fields[1]
    tx_ids = fields[3].split(',')
    if len(tx_ids) > 1:
        continue
    tx_id = tx_ids[0]
    if read_id in read_id_to_length:
        if tx_id in tx_id_to_lengths: #No need to initialize
            tx_id_to_lengths[tx_id][0] += 1
            tx_id_to_lengths[tx_id][1].append(read_id_to_length[read_id])
            gene_id_to_lengths[gene_id][0] += 1
            gene_id_to_lengths[gene_id].append(read_id_to_length[read_id])
        elif gene_id in gene_id_to_lengths: #other isoform seen?
            tx_id_to_lengths[tx_id] = [1, [read_id_to_length[read_id]]]
            gene_id_to_lengths[gene_id][0] += 1
            gene_id_to_lengths[gene_id].append(read_id_to_length[read_id])
        else: #initialize both
            tx_id_to_lengths[tx_id] = [1, [read_id_to_length[read_id]]]
            gene_id_to_lengths[gene_id] = [1, [read_id_to_length[read_id]]]

tmp_tx_id_file = open("../../results/correctionLogs/l3.gene.txt")
for line in tmp_tx_id_file:
    fields = line.strip().split()
    if len(fields) != 4:
        continue
    read_id = fields[0]
    gene_id = fields[1]
    tx_ids = fields[3].split(',')
    if len(tx_ids) > 1:
        continue
    tx_id = tx_ids[0]
    if read_id in read_id_to_length:
        if tx_id in tx_id_to_lengths: #No need to initialize
            tx_id_to_lengths[tx_id][0] += 1
            tx_id_to_lengths[tx_id][1].append(read_id_to_length[read_id])
            gene_id_to_lengths[gene_id][0] += 1
            gene_id_to_lengths[gene_id].append(read_id_to_length[read_id])
        elif gene_id in gene_id_to_lengths: #other isoform seen?
            tx_id_to_lengths[tx_id] = [1, [read_id_to_length[read_id]]]
            gene_id_to_lengths[gene_id][0] += 1
            gene_id_to_lengths[gene_id].append(read_id_to_length[read_id])
        else: #initialize both
            tx_id_to_lengths[tx_id] = [1, [read_id_to_length[read_id]]]
            gene_id_to_lengths[gene_id] = [1, [read_id_to_length[read_id]]]
            
tmp_tx_id_file = open("../../results/correctionLogs/l4.gene.txt")
for line in tmp_tx_id_file:
    fields = line.strip().split()
    if len(fields) != 4:
        continue
    read_id = fields[0]
    gene_id = fields[1]
    tx_ids = fields[3].split(',')
    if len(tx_ids) > 1:
        continue
    tx_id = tx_ids[0]
    if read_id in read_id_to_length:
        if tx_id in tx_id_to_lengths: #No need to initialize
            tx_id_to_lengths[tx_id][0] += 1
            tx_id_to_lengths[tx_id][1].append(read_id_to_length[read_id])
            gene_id_to_lengths[gene_id][0] += 1
            gene_id_to_lengths[gene_id].append(read_id_to_length[read_id])
        elif gene_id in gene_id_to_lengths: #other isoform seen?
            tx_id_to_lengths[tx_id] = [1, [read_id_to_length[read_id]]]
            gene_id_to_lengths[gene_id][0] += 1
            gene_id_to_lengths[gene_id].append(read_id_to_length[read_id])
        else: #initialize both
            tx_id_to_lengths[tx_id] = [1, [read_id_to_length[read_id]]]
            gene_id_to_lengths[gene_id] = [1, [read_id_to_length[read_id]]]
            
tmp_tx_id_file = open("../../results/correctionLogs/ya.gene.txt")
for line in tmp_tx_id_file:
    fields = line.strip().split()
    if len(fields) != 4:
        continue
    read_id = fields[0]
    gene_id = fields[1]
    tx_ids = fields[3].split(',')
    if len(tx_ids) > 1:
        continue
    tx_id = tx_ids[0]
    if read_id in read_id_to_length:
        if tx_id in tx_id_to_lengths: #No need to initialize
            tx_id_to_lengths[tx_id][0] += 1
            tx_id_to_lengths[tx_id][1].append(read_id_to_length[read_id])
            gene_id_to_lengths[gene_id][0] += 1
            gene_id_to_lengths[gene_id].append(read_id_to_length[read_id])
        elif gene_id in gene_id_to_lengths: #other isoform seen?
            tx_id_to_lengths[tx_id] = [1, [read_id_to_length[read_id]]]
            gene_id_to_lengths[gene_id][0] += 1
            gene_id_to_lengths[gene_id].append(read_id_to_length[read_id])
        else: #initialize both
            tx_id_to_lengths[tx_id] = [1, [read_id_to_length[read_id]]]
            gene_id_to_lengths[gene_id] = [1, [read_id_to_length[read_id]]]
            
tmp_tx_id_file = open("../../results/correctionLogs/ga.gene.txt")
for line in tmp_tx_id_file:
    fields = line.strip().split()
    if len(fields) != 4:
        continue
    read_id = fields[0]
    gene_id = fields[1]
    tx_ids = fields[3].split(',')
    if len(tx_ids) > 1:
        continue
    tx_id = tx_ids[0]
    if read_id in read_id_to_length:
        if tx_id in tx_id_to_lengths: #No need to initialize
            tx_id_to_lengths[tx_id][0] += 1
            tx_id_to_lengths[tx_id][1].append(read_id_to_length[read_id])
            gene_id_to_lengths[gene_id][0] += 1
            gene_id_to_lengths[gene_id].append(read_id_to_length[read_id])
        elif gene_id in gene_id_to_lengths: #other isoform seen?
            tx_id_to_lengths[tx_id] = [1, [read_id_to_length[read_id]]]
            gene_id_to_lengths[gene_id][0] += 1
            gene_id_to_lengths[gene_id].append(read_id_to_length[read_id])
        else: #initialize both
            tx_id_to_lengths[tx_id] = [1, [read_id_to_length[read_id]]]
            gene_id_to_lengths[gene_id] = [1, [read_id_to_length[read_id]]]
            
            
tmp_tx_id_file = open("../../results/correctionLogs/ml.gene.txt")
for line in tmp_tx_id_file:
    fields = line.strip().split()
    if len(fields) != 4:
        continue
    read_id = fields[0]
    gene_id = fields[1]
    tx_ids = fields[3].split(',')
    if len(tx_ids) > 1:
        continue
    tx_id = tx_ids[0]
    if read_id in read_id_to_length:
        if tx_id in tx_id_to_lengths: #No need to initialize
            tx_id_to_lengths[tx_id][0] += 1
            tx_id_to_lengths[tx_id][1].append(read_id_to_length[read_id])
            gene_id_to_lengths[gene_id][0] += 1
            gene_id_to_lengths[gene_id].append(read_id_to_length[read_id])
        elif gene_id in gene_id_to_lengths: #other isoform seen?
            tx_id_to_lengths[tx_id] = [1, [read_id_to_length[read_id]]]
            gene_id_to_lengths[gene_id][0] += 1
            gene_id_to_lengths[gene_id].append(read_id_to_length[read_id])
        else: #initialize both
            tx_id_to_lengths[tx_id] = [1, [read_id_to_length[read_id]]]
            gene_id_to_lengths[gene_id] = [1, [read_id_to_length[read_id]]]
            
tx_order = []
gene_order = []

for tx, info in tx_id_to_lengths.iteritems():
    tx_order.append((info[0],tx))


totaltxs = len(tx_order)
tx_order.sort()

for i, pair in enumerate(tx_order):
    if i > totaltxs / 3 and i <= (totaltxs * 2) / 3:
        bin_num = "Mid_expression"
    elif i > (totaltxs * 2) / 3:
        bin_num = "High_expression"
    else:
        bin_num = "Low_expression"
    tx_id = pair[1]
    median_polya = numpy.median(tx_id_to_lengths[tx_id][1])
    std_polya = numpy.std(tx_id_to_lengths[tx_id][1])

for pair in tx_order:
    tx_id = pair[1]
    tmpStr = ""
    for d in tx_id_to_lengths[tx_id][1]:
        tmpStr = "%s\t%.4f" %(tmpStr,d)
    txoutfile2.write("%s%s\n" %(tx_id,tmpStr))
