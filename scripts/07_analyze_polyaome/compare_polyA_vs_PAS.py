#!/usr/bin/env python2

import sys
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.stats

def load_read_id_to_length():
    read_id_to_length = {}### For each read in the polyA file assign its polyA tail length in a dictionary
    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L1/ce11_gen_L1_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L1/ce11_gen_L1_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L2/ce11_gen_L2_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L2/ce11_gen_L2_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])            

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L3/ce11_gen_L3_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L3/ce11_gen_L3_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L4/ce11_gen_L4_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L4/ce11_gen_L4_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])
            
    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/YA/ce11_gen_young_adult_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/YA/ce11_gen_young_adult_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/GA/ce11_gen_adult_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/GA/ce11_gen_adult_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])
            
    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/ML/ce11_gen_male_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/ML/ce11_gen_male_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])
    return read_id_to_length

cluster_assignment_file = open(sys.argv[1],'r')
pas_assignment_file = open(sys.argv[2],'r')
pas_to_polyA_out = open(sys.argv[3],'w')

read_id_to_length = load_read_id_to_length()
cluster_to_lengths = {}
cluster_to_PAS = {}
PAS_to_lengths = {}

for line in cluster_assignment_file:
    fields = line.strip().split()
    read_id = fields[0]
    cluster_id = fields[1]
    if read_id in read_id_to_length:
        length = read_id_to_length[read_id]
    else:
        continue
    if cluster_id in cluster_to_lengths:
        cluster_to_lengths[cluster_id].append(length)
    else:
        cluster_to_lengths[cluster_id] = [length]
    
for line in pas_assignment_file:
    fields = line.strip().split()
    cluster_id = fields[0]
    pas = fields[1]
    if pas != "AATAAA" and pas != "noPAS":
    #if pas != "AATAAA" and pas != "noPAS" and pas != "AATGAA" and pas != "TATAAA" and pas != "CATAAA":
    #if False:
        pas_group = "altPAS"
        #pas_group = "otherAltPAS"
    else:
        pas_group = pas
    cluster_to_PAS[cluster_id] = pas_group

for cluster_id in cluster_to_PAS:
    pas_group = cluster_to_PAS[cluster_id]
    if cluster_id in cluster_to_lengths:
        lengths = cluster_to_lengths[cluster_id]
        if pas_group in PAS_to_lengths:
            for length in lengths:
                PAS_to_lengths[pas_group].append(length)
        else:
            PAS_to_lengths[pas_group] = lengths
pas_to_polyA_out.write("length\tPAS\n")

# print "canon vs alt"
# print scipy.stats.mannwhitneyu(PAS_to_lengths["AATAAA"],PAS_to_lengths["altPAS"])
# print "canon vs no"
# print scipy.stats.mannwhitneyu(PAS_to_lengths["AATAAA"],PAS_to_lengths["noPAS"])
# print "alt vs no"
# print scipy.stats.mannwhitneyu(PAS_to_lengths["altPAS"],PAS_to_lengths["noPAS"])


for pas_group in PAS_to_lengths:
    for length in PAS_to_lengths[pas_group]:
        pas_to_polyA_out.write("%f\t%s\n" %(length,pas_group))
