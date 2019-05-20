#!/usr/bin/env python2

import sets
import pandas
import numpy
import seaborn
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hc
import scipy.spatial as sp


seaborn.set(font_scale=2.5)

read_to_cluster_in = open("../../results/utrs/assignments/all_isoforms_utrs.tsv",'r')
bed_in = open("../../results/utrs/beds/all_isoforms_utrs.bed",'r')

l1_in = open("../../results/isoforms/L1_isoforms.tsv",'r')
l2_in = open("../../results/isoforms/L2_isoforms.tsv",'r')
l3_in = open("../../results/isoforms/L3_isoforms.tsv",'r')
l4_in = open("../../results/isoforms/L4_isoforms.tsv",'r')
ya_in = open("../../results/isoforms/YA_isoforms.tsv",'r')
ga_in = open("../../results/isoforms/GA_isoforms.tsv",'r')
ml_in = open("../../results/isoforms/ML_isoforms.tsv",'r')

l1_bed_out = open("../../results/utrs/beds/L1_utrs.bed",'w')
l2_bed_out = open("../../results/utrs/beds/L2_utrs.bed",'w')
l3_bed_out = open("../../results/utrs/beds/L3_utrs.bed",'w')
l4_bed_out = open("../../results/utrs/beds/L4_utrs.bed",'w')
ya_bed_out = open("../../results/utrs/beds/YA_utrs.bed",'w')
ga_bed_out = open("../../results/utrs/beds/GA_utrs.bed",'w')
ml_bed_out = open("../../results/utrs/beds/ML_utrs.bed",'w')

l1_tsv_out = open("../../results/utrs/assignments/L1_utrs.tsv",'w')
l2_tsv_out = open("../../results/utrs/assignments/L2_utrs.tsv",'w')
l3_tsv_out = open("../../results/utrs/assignments/L3_utrs.tsv",'w')
l4_tsv_out = open("../../results/utrs/assignments/L4_utrs.tsv",'w')
ya_tsv_out = open("../../results/utrs/assignments/YA_utrs.tsv",'w')
ga_tsv_out = open("../../results/utrs/assignments/GA_utrs.tsv",'w')
ml_tsv_out = open("../../results/utrs/assignments/ML_utrs.tsv",'w')

l1_clusters = sets.Set()
l2_clusters = sets.Set()
l3_clusters = sets.Set()
l4_clusters = sets.Set()
ya_clusters = sets.Set()
ga_clusters = sets.Set()
ml_clusters = sets.Set()
total_clusters  = sets.Set()

read_to_cluster_dict = {}

for line in read_to_cluster_in:
    fields = line.strip().split()
    read_id = fields[0]
    cluster_id = fields[1]
    read_to_cluster_dict[read_id] = cluster_id
    total_clusters.add(cluster_id)

print len(total_clusters)

bed_dict = {}
for line in bed_in:
    fields = line.strip().split()
    cluster_id = fields[3]
    bed_dict[cluster_id] = line

#l1
for line in l1_in:
    fields = line.strip().split()
    read_ids = fields[13].split(',')
    for read_id in read_ids:
        if read_id in read_to_cluster_dict: #if read has an assigned utr cluster:
            cluster_id = read_to_cluster_dict[read_id] #fetch assigned cluster
            l1_tsv_out.write("%s\t%s\n" %(read_id,cluster_id)) #write read & cluster assignment to stage specific assignment file
            l1_clusters.add(cluster_id) # add cluster to set of observed l1 clusters

for cluster_id in l1_clusters:
    l1_bed_out.write(bed_dict[cluster_id]) #write cluster to stage specific 3'UTR bed file 

#l2
for line in l2_in:
    fields = line.strip().split()
    read_ids = fields[13].split(',')
    for read_id in read_ids:
        if read_id in read_to_cluster_dict: #if read has an assigned utr cluster:
            cluster_id = read_to_cluster_dict[read_id] #fetch assigned cluster
            l2_tsv_out.write("%s\t%s\n" %(read_id,cluster_id)) #write read & cluster assignment to stage specific assignment file
            l2_clusters.add(cluster_id) # add cluster to set of observed l2 clusters

for cluster_id in l2_clusters:
    l2_bed_out.write(bed_dict[cluster_id]) #write cluster to stage specific 3'UTR bed file 

#l3
for line in l3_in:
    fields = line.strip().split()
    read_ids = fields[13].split(',')
    for read_id in read_ids:
        if read_id in read_to_cluster_dict: #if read has an assigned utr cluster:
            cluster_id = read_to_cluster_dict[read_id] #fetch assigned cluster
            l3_tsv_out.write("%s\t%s\n" %(read_id,cluster_id)) #write read & cluster assignment to stage specific assignment file
            l3_clusters.add(cluster_id) # add cluster to set of observed l3 clusters

for cluster_id in l3_clusters:
    l3_bed_out.write(bed_dict[cluster_id]) #write cluster to stage specific 3'UTR bed file 

#l4
for line in l4_in:
    fields = line.strip().split()
    read_ids = fields[13].split(',')
    for read_id in read_ids:
        if read_id in read_to_cluster_dict: #if read has an assigned utr cluster:
            cluster_id = read_to_cluster_dict[read_id] #fetch assigned cluster
            l4_tsv_out.write("%s\t%s\n" %(read_id,cluster_id)) #write read & cluster assignment to stage specific assignment file
            l4_clusters.add(cluster_id) # add cluster to set of observed l4 clusters

for cluster_id in l4_clusters:
    l4_bed_out.write(bed_dict[cluster_id]) #write cluster to stage specific 3'UTR bed file 

#ya
for line in ya_in:
    fields = line.strip().split()
    read_ids = fields[13].split(',')
    for read_id in read_ids:
        if read_id in read_to_cluster_dict: #if read has an assigned utr cluster:
            cluster_id = read_to_cluster_dict[read_id] #fetch assigned cluster
            ya_tsv_out.write("%s\t%s\n" %(read_id,cluster_id)) #write read & cluster assignment to stage specific assignment file
            ya_clusters.add(cluster_id) # add cluster to set of observed ya clusters

for cluster_id in ya_clusters:
    ya_bed_out.write(bed_dict[cluster_id]) #write cluster to stage specific 3'UTR bed file 

#ga
for line in ga_in:
    fields = line.strip().split()
    read_ids = fields[13].split(',')
    for read_id in read_ids:
        if read_id in read_to_cluster_dict: #if read has an assigned utr cluster:
            cluster_id = read_to_cluster_dict[read_id] #fetch assigned cluster
            ga_tsv_out.write("%s\t%s\n" %(read_id,cluster_id)) #write read & cluster assignment to stage specific assignment file
            ga_clusters.add(cluster_id) # add cluster to set of observed ga clusters

for cluster_id in ga_clusters:
    ga_bed_out.write(bed_dict[cluster_id]) #write cluster to stage specific 3'UTR bed file 

#ml
for line in ml_in:
    fields = line.strip().split()
    read_ids = fields[13].split(',')
    for read_id in read_ids:
        if read_id in read_to_cluster_dict: #if read has an assigned utr cluster:
            cluster_id = read_to_cluster_dict[read_id] #fetch assigned cluster
            ml_tsv_out.write("%s\t%s\n" %(read_id,cluster_id)) #write read & cluster assignment to stage specific assignment file
            ml_clusters.add(cluster_id) # add cluster to set of observed ml clusters

for cluster_id in ml_clusters:
    ml_bed_out.write(bed_dict[cluster_id]) #write cluster to stage specific 3'UTR bed file 




utr_outfile = open("../../results/overlaps/overlap.utr.matrix",'w') 

all_clusters = l1_clusters.union(l2_clusters.union(l3_clusters.union(l4_clusters.union(ya_clusters.union(ga_clusters.union(ml_clusters))))))

l1_uniq = open("../../results/utrs/beds/uniqs/L1_uniq.bed",'w')
l2_uniq = open("../../results/utrs/beds/uniqs/L2_uniq.bed",'w')
l3_uniq = open("../../results/utrs/beds/uniqs/L3_uniq.bed",'w')
l4_uniq = open("../../results/utrs/beds/uniqs/L4_uniq.bed",'w')
ya_uniq = open("../../results/utrs/beds/uniqs/YA_uniq.bed",'w')
ga_uniq = open("../../results/utrs/beds/uniqs/GA_uniq.bed",'w')
ml_uniq = open("../../results/utrs/beds/uniqs/ML_uniq.bed",'w')
for cluster_id in all_clusters:
    utr_outfile.write("%i\t%i\t%i\t%i\t%i\t%i\t%i\n" %(cluster_id in l1_clusters,cluster_id in l2_clusters,cluster_id in l3_clusters,cluster_id in l4_clusters,cluster_id in ya_clusters,cluster_id in ga_clusters,cluster_id in ml_clusters))
    ## Grab all the unique UTRs
    if cluster_id in l1_clusters and cluster_id not in l2_clusters and cluster_id not in l3_clusters and cluster_id not in l4_clusters and cluster_id not in ya_clusters and cluster_id not in ga_clusters and cluster_id not in ml_clusters:
        l1_uniq.write(bed_dict[cluster_id])
    if cluster_id not in l1_clusters and cluster_id in l2_clusters and cluster_id not in l3_clusters and cluster_id not in l4_clusters and cluster_id not in ya_clusters and cluster_id not in ga_clusters and cluster_id not in ml_clusters:
        l2_uniq.write(bed_dict[cluster_id])
    if cluster_id not in l1_clusters and cluster_id not in l2_clusters and cluster_id in l3_clusters and cluster_id not in l4_clusters and cluster_id not in ya_clusters and cluster_id not in ga_clusters and cluster_id not in ml_clusters:
        l3_uniq.write(bed_dict[cluster_id])
    if cluster_id not in l1_clusters and cluster_id not in l2_clusters and cluster_id not in l3_clusters and cluster_id in l4_clusters and cluster_id not in ya_clusters and cluster_id not in ga_clusters and cluster_id not in ml_clusters:
        l4_uniq.write(bed_dict[cluster_id])
    if cluster_id not in l1_clusters and cluster_id not in l2_clusters and cluster_id not in l3_clusters and cluster_id not in l4_clusters and cluster_id in ya_clusters and cluster_id not in ga_clusters and cluster_id not in ml_clusters:
        ya_uniq.write(bed_dict[cluster_id])
    if cluster_id not in l1_clusters and cluster_id not in l2_clusters and cluster_id not in l3_clusters and cluster_id not in l4_clusters and cluster_id not in ya_clusters and cluster_id in ga_clusters and cluster_id not in ml_clusters:
        ga_uniq.write(bed_dict[cluster_id])
    if cluster_id not in l1_clusters and cluster_id not in l2_clusters and cluster_id not in l3_clusters and cluster_id not in l4_clusters and cluster_id not in ya_clusters and cluster_id not in ga_clusters and cluster_id in ml_clusters:
        ml_uniq.write(bed_dict[cluster_id])

#Plot heatmap
# set_list = [l1_clusters,l2_clusters,l3_clusters,l4_clusters,ya_clusters,ga_clusters,ml_clusters]
# stages = ["L1","L2","L3","L4","YA","GA","ML"]
# heat = pandas.DataFrame(numpy.ones([7,7]),columns=["L1","L2","L3","L4","young adult","mature adult","male"])
# heat["labels"] = ["L1","L2","L3","L4","young adult","mature adult","male"]
# heat = heat.set_index("labels")
# labels = ["L1","L2","L3","L4","young adult","mature adult","male"]
# set_list = [l1_clusters,l2_clusters,l3_clusters,l4_clusters,ya_clusters,ga_clusters]
# stages = ["L1","L2","L3","L4","YA","GA","ML"]
# heat = pandas.DataFrame(numpy.ones([6,6]),columns=["L1","L2","L3","L4","young adult","mature adult",])
# heat["labels"] = ["L1","L2","L3","L4","young adult","mature adult",]
# heat = heat.set_index("labels")
# labels = ["L1","L2","L3","L4","young adult","mature adult"]
# for x in range(len(set_list)):
#     for y in range(len(set_list)):
#         lab1 = labels[x]
#         lab2 = labels[y]
#         temp = float(len(set_list[x].intersection(set_list[y]))) / float(len(set_list[x].union(set_list[y])))
#         heat[lab1][lab2] = temp
#         heat[lab2][lab1] = temp
#
# row_linkage,col_linkage = (hc.linkage(sp.distance.pdist(x),method="average") for x in (heat.values,heat.values.T))
#
# row_linkage2 = row_linkage
#
# tmp = row_linkage2[2][0]
# row_linkage2[2][0] = row_linkage2[2][1]
# row_linkage2[2][1] = tmp
# #
# tmp = row_linkage2[3][0]
# row_linkage2[3][0] = row_linkage2[3][1]
# row_linkage2[3][1] = tmp
# #
# tmp = row_linkage2[4][0]
# row_linkage2[4][0] = row_linkage2[4][1]
# row_linkage2[4][1] = tmp
#
# col_linkage2 = col_linkage
# #
# tmp = col_linkage2[2][0]
# col_linkage2[2][0] = col_linkage2[2][1]
# col_linkage2[2][1] = tmp
# #
# tmp = col_linkage2[3][0]
# col_linkage2[3][0] = col_linkage2[3][1]
# col_linkage2[3][1] = tmp
# #
# tmp = col_linkage2[4][0]
# col_linkage2[4][0] = col_linkage2[4][1]
# col_linkage2[4][1] = tmp
#
# seaborn.clustermap(heat,vmin=0.4,vmax=1.0,cmap="viridis",row_linkage=row_linkage2,col_linkage=col_linkage2)
#
# #plt.show()
# #seaborn.heatmap(heat,cmap="viridis")
# plt.title("UTR Jaccard Indices")
# plt.savefig("plots/heatmaps/lowercase_utr_heatmap.pdf",bbox_inches="tight")
# plt.clf()

print "l1 - %d utrs" %(len(l1_clusters))
print "l2 - %d utrs" %(len(l2_clusters))
print "l3 - %d utrs" %(len(l3_clusters))
print "l4 - %d utrs" %(len(l4_clusters))
print "ya - %d utrs" %(len(ya_clusters))
print "ga - %d utrs" %(len(ga_clusters))
print "ml - %d utrs" %(len(ml_clusters))
print "total - %d utrs" %(len(all_clusters))

