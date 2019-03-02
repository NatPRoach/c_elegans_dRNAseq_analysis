#!/usr/bin/env python2

import sets

novel_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/novel/all_novel_utrs.bed")
l1_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/L1_utrs.bed")
l2_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/L2_utrs.bed")
l3_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/L3_utrs.bed")
l4_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/L4_utrs.bed")
ya_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/YA_utrs.bed")
ga_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/GA_utrs.bed")
ml_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/ML_utrs.bed")

l1_out = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/novel/L1_novel_utrs.bed",'w')
l2_out = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/novel/L2_novel_utrs.bed",'w')
l3_out = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/novel/L3_novel_utrs.bed",'w')
l4_out = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/novel/L4_novel_utrs.bed",'w')
ya_out = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/novel/YA_novel_utrs.bed",'w')
ga_out = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/novel/GA_novel_utrs.bed",'w')
ml_out = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/novel/ML_novel_utrs.bed",'w')

outfile = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/analyzeUTRome/novelUTRs/novelUTRs.txt",'w')

novel_set = sets.Set()
l1_set = sets.Set()
l2_set = sets.Set()
l3_set = sets.Set()
l4_set = sets.Set()
ya_set = sets.Set()
ga_set = sets.Set()
ml_set = sets.Set()

novel_gene_set = sets.Set()
l1_gene_set = sets.Set()
l2_gene_set = sets.Set()
l3_gene_set = sets.Set()
l4_gene_set = sets.Set()
ya_gene_set = sets.Set()
ga_gene_set = sets.Set()
ml_gene_set = sets.Set()

cluster_to_line = {}
for line in novel_in:
    fields = line.strip().split()
    cluster_id = fields[3]
    novel_set.add(cluster_id)
    gene_id = fields[3].split('-')[0]
    novel_gene_set.add(gene_id)
    cluster_to_line[cluster_id] = line
for line in l1_in:
    fields = line.strip().split()
    cluster_id = fields[3]
    l1_set.add(cluster_id)
    # gene_id = fields[3].split('-')[0]
    # l1_gene_set.add(gene_id)
for line in l2_in:
    fields = line.strip().split()
    cluster_id = fields[3]
    l2_set.add(cluster_id)
    # gene_id = fields[3].split('-')[0]
    # l2_gene_set.add(gene_id)
for line in l3_in:
    fields = line.strip().split()
    cluster_id = fields[3]
    l3_set.add(cluster_id)
    # gene_id = fields[3].split('-')[0]
    # l3_gene_set.add(gene_id)
for line in l4_in:
    fields = line.strip().split()
    cluster_id = fields[3]
    l4_set.add(cluster_id)
    # gene_id = fields[3].split('-')[0]
    # l4_gene_set.add(gene_id)
for line in ya_in:
    fields = line.strip().split()
    cluster_id = fields[3]
    ya_set.add(cluster_id)
    # gene_id = fields[3].split('-')[0]
    # ya_gene_set.add(gene_id)
for line in ga_in:
    fields = line.strip().split()
    cluster_id = fields[3]
    ga_set.add(cluster_id)
    # gene_id = fields[3].split('-')[0]
    # ga_gene_set.add(gene_id)
for line in ml_in:
    fields = line.strip().split()
    cluster_id = fields[3]
    ml_set.add(cluster_id)
    # gene_id = fields[3].split('-')[0]
    # ml_gene_set.add(gene_id)

l1_novel = novel_set.intersection(l1_set)
l2_novel = novel_set.intersection(l2_set)
l3_novel = novel_set.intersection(l3_set)
l4_novel = novel_set.intersection(l4_set)
ya_novel = novel_set.intersection(ya_set)
ga_novel = novel_set.intersection(ga_set)
ml_novel = novel_set.intersection(ml_set)

for cluster in l1_novel:
    l1_out.write(cluster_to_line[cluster])
    gene_id = cluster.split('-')[0]
    l1_gene_set.add(gene_id)
for cluster in l2_novel:
    l2_out.write(cluster_to_line[cluster])
    gene_id = cluster.split('-')[0]
    l2_gene_set.add(gene_id)
for cluster in l3_novel:
    l3_out.write(cluster_to_line[cluster])
    gene_id = cluster.split('-')[0]
    l3_gene_set.add(gene_id)
for cluster in l4_novel:
    l4_out.write(cluster_to_line[cluster])
    gene_id = cluster.split('-')[0]
    l4_gene_set.add(gene_id)
for cluster in ya_novel:
    ya_out.write(cluster_to_line[cluster])
    gene_id = cluster.split('-')[0]
    ya_gene_set.add(gene_id)
for cluster in ga_novel:
    ga_out.write(cluster_to_line[cluster])
    gene_id = cluster.split('-')[0]
    ga_gene_set.add(gene_id)
for cluster in ml_novel:
    ml_out.write(cluster_to_line[cluster])
    gene_id = cluster.split('-')[0]
    ml_gene_set.add(gene_id)
    
#
# print len(novel_set)
# print len(l1_novel)
# print len(l2_novel)
# print len(l3_novel)
# print len(l4_novel)
# print len(ya_novel)
# print len(ga_novel)
# print len(ml_novel)
# print ""
# print len(novel_gene_set)
# print len(l1_gene_set)
# print len(l2_gene_set)
# print len(l3_gene_set)
# print len(l4_gene_set)
# print len(ya_gene_set)
# print len(ga_gene_set)
# print len(ml_gene_set)

outfile.write("stage\tdataset\tcounts\tx_order\n")

outfile.write("L1\tnovel UTRs\t%d\t0\n" %(len(l1_novel)))
outfile.write("L2\tnovel UTRs\t%d\t1\n" %(len(l2_novel)))
outfile.write("L3\tnovel UTRs\t%d\t2\n" %(len(l3_novel)))
outfile.write("L4\tnovel UTRs\t%d\t3\n" %(len(l4_novel)))
outfile.write("young adult\tnovel UTRs\t%d\t4\n" %(len(ya_novel)))
outfile.write("mature adult\tnovel UTRs\t%d\t5\n" %(len(ga_novel)))
outfile.write("male\tnovel UTRs\t%d\t6\n" %(len(ml_novel)))
outfile.write("all\tnovel UTRs\t%d\t7\n" %(len(novel_set)))

outfile.write("L1\tgenes with novel UTRs\t%d\t0\n" %(len(l1_gene_set)))
outfile.write("L2\tgenes with novel UTRs\t%d\t1\n" %(len(l2_gene_set)))
outfile.write("L3\tgenes with novel UTRs\t%d\t2\n" %(len(l3_gene_set)))
outfile.write("L4\tgenes with novel UTRs\t%d\t3\n" %(len(l4_gene_set)))
outfile.write("young adult\tgenes with novel UTRs\t%d\t4\n" %(len(ya_gene_set)))
outfile.write("mature adult\tgenes with novel UTRs\t%d\t5\n" %(len(ga_gene_set)))
outfile.write("male\tgenes with novel UTRs\t%d\t6\n" %(len(ml_gene_set)))
outfile.write("all\tgenes with novel UTRs\t%d\t7\n" %(len(novel_gene_set)))