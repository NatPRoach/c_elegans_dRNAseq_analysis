#!/usr/bin/env python
import sys
import sets

def load_read_id_to_length():
    read_id_to_length = {}### For each read in the polyA file assign its polyA tail length in a dictionary
    tmppolyafile = open("../../data/L1/bio1/tech1/analysis/ce11_gen_L1_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0].upper()] = float(fields[8])

    tmppolyafile = open("../../data/L1/bio1/tech2/analysis/ce11_gen_L1_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0].upper()] = float(fields[8])

    tmppolyafile = open("../../data/L2/bio1/tech1/analysis/ce11_gen_L2_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0].upper()] = float(fields[8])

    tmppolyafile = open("../../data/L2/bio1/tech2/analysis/ce11_gen_L2_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0].upper()] = float(fields[8])            

    tmppolyafile = open("../../data/L3/bio1/tech1/analysis/ce11_gen_L3_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0].upper()] = float(fields[8])

    tmppolyafile = open("../../data/L3/bio1/tech2/analysis/ce11_gen_L3_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0].upper()] = float(fields[8])

    tmppolyafile = open("../../data/L4/bio1/tech1/analysis/ce11_gen_L4_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0].upper()] = float(fields[8])

    tmppolyafile = open("../../data/L4/bio1/tech2/analysis/ce11_gen_L4_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0].upper()] = float(fields[8])
            
    tmppolyafile = open("../../data/young_adult/bio1/tech1/analysis/ce11_gen_young_adult_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0].upper()] = float(fields[8])

    tmppolyafile = open("../../data/young_adult/bio1/tech2/analysis/ce11_gen_young_adult_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0].upper()] = float(fields[8])

    tmppolyafile = open("../../data/adult/bio1/tech1/analysis/ce11_gen_adult_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0].upper()] = float(fields[8])

    tmppolyafile = open("../../data/adult/bio1/tech2/analysis/ce11_gen_adult_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0].upper()] = float(fields[8])
            
    tmppolyafile = open("../../data/male/bio1/tech1/analysis/ce11_gen_male_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0].upper()] = float(fields[8])

    tmppolyafile = open("../../data/male/bio1/tech2/analysis/ce11_gen_male_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0].upper()] = float(fields[8])
    return read_id_to_length


# infile = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/isoforms/all_stringent_isoforms.tsv")
infile0 = open(sys.argv[1])
count = 0
non_coding_reads = sets.Set()
coding_reads = sets.Set()
for i,line in enumerate(infile0):
    if i == 0:
        continue
    fields = line.strip().split()
    if float(fields[-1]) < 0.5:
        gene,read_ids = fields[0].split('_')
        for read_id in read_ids.split(','):
            non_coding_reads.add(read_id.upper())
        count += 1
    else:
        gene,read_ids = fields[0].split('_')
        for read_id in read_ids.split(','):
            coding_reads.add(read_id.upper())




# infile1=open("/Users/nproach/Documents/c_elegans_dRNAseq_analysis/results/scratch/PASanalysis/assignmentFiles/L1_stringent.PASassignments.txt")
infile1 = open(sys.argv[2])
beta_canon_clusters = sets.Set()
beta_alt_clusters = sets.Set()
beta_no_clusters = sets.Set()

for line in infile1:
    fields = line.strip().split()
    if fields[1] == 'AATAAA':
        beta_canon_clusters.add(fields[0])
    elif fields[1] == 'noPAS':
        beta_no_clusters.add(fields[0])
    else:
        beta_alt_clusters.add(fields[0])

beta_canon_reads = sets.Set()
beta_alt_reads = sets.Set()
beta_no_reads = sets.Set()



# infile2=open("/Users/nproach/Documents/c_elegans_dRNAseq_analysis/results/utrs/assignments/L1_stringent_utrs.tsv")
infile2 = open(sys.argv[3])
for line in infile2:
    fields = line.strip().split()
    if fields[1] in beta_canon_clusters:
        beta_canon_reads.add(fields[0].upper())
    elif fields[1] in beta_alt_clusters:
        beta_alt_reads.add(fields[0].upper())
    elif fields[1] in beta_no_clusters:
        beta_no_reads.add(fields[0].upper())

read_id_to_length = load_read_id_to_length()

outfile = open(sys.argv[4],'w')
outfile.write("length\tCoding\tPAS\n")
for read in beta_canon_reads:
    if read in read_id_to_length:
        length = read_id_to_length[read]
        if read in coding_reads:
            coding = "Coding"
        elif read in non_coding_reads:
            coding = "Non-coding"
        else:
            continue
        # print "%f\t%s\t%s" %(length,coding,"AAUAAA")
        outfile.write("%f\t%s\t%s\n" %(length,coding,"AATAAA"))
        

for read in beta_alt_reads:
    if read in read_id_to_length:
        length = read_id_to_length[read]
        if read in coding_reads:
            coding = "Coding"
        elif read in non_coding_reads:
            coding = "Non-coding"
        else:
            continue
        # print "%f\t%s\t%s" %(length,coding,"AltPAS")
        outfile.write("%f\t%s\t%s\n" %(length,coding,"altPAS"))
        

for read in beta_no_reads:
    if read in read_id_to_length:
        length = read_id_to_length[read]
        if read in coding_reads:
            coding = "Coding"
        elif read in non_coding_reads:
            coding = "Non-coding"
        else:
            continue
        # print "%f\t%s\t%s" %(length,coding,"NoPAS")
        outfile.write("%f\t%s\t%s\n" %(length,coding,"noPAS"))
        



# import sys
# import matplotlib.pyplot as plt
# import seaborn as sns
# import scipy.stats
#
# def load_read_id_to_length():
#     read_id_to_length = {}### For each read in the polyA file assign its polyA tail length in a dictionary
#     tmppolyafile = open("../../data/L1/bio1/tech1/analysis/ce11_gen_L1_bio1_tech1_Minimap2_k14.polya",'r')
#     for i, line in enumerate(tmppolyafile):
#         if i != 0:
#             fields = line.strip().split()
#             if fields[9] == "PASS":
#                 read_id_to_length[fields[0]] = float(fields[8])
#
#     tmppolyafile = open("../../data/L1/bio1/tech2/analysis/ce11_gen_L1_bio1_tech2_Minimap2_k14.polya",'r')
#     for i, line in enumerate(tmppolyafile):
#         if i != 0:
#             fields = line.strip().split()
#             if fields[9] == "PASS":
#                 read_id_to_length[fields[0]] = float(fields[8])
#
#     tmppolyafile = open("../../data/L2/bio1/tech1/analysis/ce11_gen_L2_bio1_tech1_Minimap2_k14.polya",'r')
#     for i, line in enumerate(tmppolyafile):
#         if i != 0:
#             fields = line.strip().split()
#             if fields[9] == "PASS":
#                 read_id_to_length[fields[0]] = float(fields[8])
#
#     tmppolyafile = open("../../data/L2/bio1/tech2/analysis/ce11_gen_L2_bio1_tech2_Minimap2_k14.polya",'r')
#     for i, line in enumerate(tmppolyafile):
#         if i != 0:
#             fields = line.strip().split()
#             if fields[9] == "PASS":
#                 read_id_to_length[fields[0]] = float(fields[8])
#
#     tmppolyafile = open("../../data/L3/bio1/tech1/analysis/ce11_gen_L3_bio1_tech1_Minimap2_k14.polya",'r')
#     for i, line in enumerate(tmppolyafile):
#         if i != 0:
#             fields = line.strip().split()
#             if fields[9] == "PASS":
#                 read_id_to_length[fields[0]] = float(fields[8])
#
#     tmppolyafile = open("../../data/L3/bio1/tech2/analysis/ce11_gen_L3_bio1_tech2_Minimap2_k14.polya",'r')
#     for i, line in enumerate(tmppolyafile):
#         if i != 0:
#             fields = line.strip().split()
#             if fields[9] == "PASS":
#                 read_id_to_length[fields[0]] = float(fields[8])
#
#     tmppolyafile = open("../../data/L4/bio1/tech1/analysis/ce11_gen_L4_bio1_tech1_Minimap2_k14.polya",'r')
#     for i, line in enumerate(tmppolyafile):
#         if i != 0:
#             fields = line.strip().split()
#             if fields[9] == "PASS":
#                 read_id_to_length[fields[0]] = float(fields[8])
#
#     tmppolyafile = open("../../data/L4/bio1/tech2/analysis/ce11_gen_L4_bio1_tech2_Minimap2_k14.polya",'r')
#     for i, line in enumerate(tmppolyafile):
#         if i != 0:
#             fields = line.strip().split()
#             if fields[9] == "PASS":
#                 read_id_to_length[fields[0]] = float(fields[8])
#
#     tmppolyafile = open("../../data/young_adult/bio1/tech1/analysis/ce11_gen_young_adult_bio1_tech1_Minimap2_k14.polya",'r')
#     for i, line in enumerate(tmppolyafile):
#         if i != 0:
#             fields = line.strip().split()
#             if fields[9] == "PASS":
#                 read_id_to_length[fields[0]] = float(fields[8])
#
#     tmppolyafile = open("../../data/young_adult/bio1/tech2/analysis/ce11_gen_young_adult_bio1_tech2_Minimap2_k14.polya",'r')
#     for i, line in enumerate(tmppolyafile):
#         if i != 0:
#             fields = line.strip().split()
#             if fields[9] == "PASS":
#                 read_id_to_length[fields[0]] = float(fields[8])
#
#     tmppolyafile = open("../../data/adult/bio1/tech1/analysis/ce11_gen_adult_bio1_tech1_Minimap2_k14.polya",'r')
#     for i, line in enumerate(tmppolyafile):
#         if i != 0:
#             fields = line.strip().split()
#             if fields[9] == "PASS":
#                 read_id_to_length[fields[0]] = float(fields[8])
#
#     tmppolyafile = open("../../data/adult/bio1/tech2/analysis/ce11_gen_adult_bio1_tech2_Minimap2_k14.polya",'r')
#     for i, line in enumerate(tmppolyafile):
#         if i != 0:
#             fields = line.strip().split()
#             if fields[9] == "PASS":
#                 read_id_to_length[fields[0]] = float(fields[8])
#
#     tmppolyafile = open("../../data/male/bio1/tech1/analysis/ce11_gen_male_bio1_tech1_Minimap2_k14.polya",'r')
#     for i, line in enumerate(tmppolyafile):
#         if i != 0:
#             fields = line.strip().split()
#             if fields[9] == "PASS":
#                 read_id_to_length[fields[0]] = float(fields[8])
#
#     tmppolyafile = open("../../data/male/bio1/tech2/analysis/ce11_gen_male_bio1_tech2_Minimap2_k14.polya",'r')
#     for i, line in enumerate(tmppolyafile):
#         if i != 0:
#             fields = line.strip().split()
#             if fields[9] == "PASS":
#                 read_id_to_length[fields[0]] = float(fields[8])
#     return read_id_to_length
#
# cluster_assignment_file = open(sys.argv[1],'r')
# pas_assignment_file = open(sys.argv[2],'r')
# pas_to_polyA_out = open(sys.argv[3],'w')
#
# read_id_to_length = load_read_id_to_length()
# cluster_to_lengths = {}
# cluster_to_PAS = {}
# PAS_to_lengths = {}
#
# for line in cluster_assignment_file:
#     fields = line.strip().split()
#     read_id = fields[0]
#     cluster_id = fields[1]
#     if read_id in read_id_to_length:
#         length = read_id_to_length[read_id]
#     else:
#         continue
#     if cluster_id in cluster_to_lengths:
#         cluster_to_lengths[cluster_id].append(length)
#     else:
#         cluster_to_lengths[cluster_id] = [length]
#
# for line in pas_assignment_file:
#     fields = line.strip().split()
#     cluster_id = fields[0]
#     pas = fields[1]
#     if pas != "AATAAA" and pas != "noPAS":
#         pas_group = "altPAS"
#     else:
#         pas_group = pas
#     cluster_to_PAS[cluster_id] = pas_group
#
# for cluster_id in cluster_to_PAS:
#     pas_group = cluster_to_PAS[cluster_id]
#     if cluster_id in cluster_to_lengths:
#         lengths = cluster_to_lengths[cluster_id]
#         if pas_group in PAS_to_lengths:
#             for length in lengths:
#                 PAS_to_lengths[pas_group].append(length)
#         else:
#             PAS_to_lengths[pas_group] = lengths
# pas_to_polyA_out.write("length\tPAS\n")
#
# # print "canon vs alt"
# # print scipy.stats.mannwhitneyu(PAS_to_lengths["AATAAA"],PAS_to_lengths["altPAS"])
# # print "canon vs no"
# # print scipy.stats.mannwhitneyu(PAS_to_lengths["AATAAA"],PAS_to_lengths["noPAS"])
# # print "alt vs no"
# # print scipy.stats.mannwhitneyu(PAS_to_lengths["altPAS"],PAS_to_lengths["noPAS"])
#
#
# for pas_group in PAS_to_lengths:
#     for length in PAS_to_lengths[pas_group]:
#         pas_to_polyA_out.write("%f\t%s\n" %(length,pas_group))
