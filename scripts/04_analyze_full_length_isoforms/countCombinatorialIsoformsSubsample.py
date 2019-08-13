#!/usr/bin/env python2

import sets
import numpy as np

### nroach2@jhu.edu
### Need to subsample the number of reads we consider in each dataset to get a sense of how many genes and isoforms we would identify if we were at saturation.
def countCombinatorialIsoforms(infilepath,utr_assignments,fraction_retained=1.0):
    read_subsample = sets.Set()
    infile = open(infilepath)
    for line in infile:
        fields = line.strip().split()
        read_ids = fields[13].split(',')
        for read_id in read_ids:
            if np.random.random() < fraction_retained:
                read_subsample.add(read_id)
    # print len(read_subsample)
    isoform_count = 0
    splice_isoform_count = 0
    total_uniq_utrs = sets.Set()
    infile = open(infilepath)
    for line in infile:
        fields = line.strip().split()
        gene = fields[0]
        chrom = fields[1]
        strand = fields[2]
        cds = fields[3]
        seq = fields[4]
        cds_start = int(fields[5])
        cds_end = int(fields[6])
        num_reads = int(fields[7])
        num_exons = int(fields[8])
        start_codon = int(fields[9])
        stop_codon = int(fields[10])
        stop_positions = [int(x) for x in fields[11].split(',')]
        read_ids = fields[13].split(',')
        read_ids2 = []
        for read_id in read_ids:
            if read_id in read_subsample:
                read_ids2.append(read_id)
        if len(read_ids2) != 0:
            splice_isoform_count += 1
        else:
            continue
        retained_intron = fields[14] == '1'
        assert len(stop_positions) == len(read_ids)
        if retained_intron: #For now ignore retained intron transcripts
            continue
        #Each line represents an intron chain (with 5' tolerance... remove 5' tolerance?)
        #Need to count number of utrs represented in each intron chain
        #We have a file of read_id -> gene & UTR call
        uniq_utrs = sets.Set()
        for read_id in read_ids2:
            if read_id in utr_assignments:
                gene_id,cluster_id = utr_assignments[read_id]
                assert gene_id == gene
                uniq_utrs.add(cluster_id)
                total_uniq_utrs.add(cluster_id)
        if len(uniq_utrs) == 0:
            utr_count = 1
        else:
            utr_count = len(uniq_utrs)
        isoform_count += utr_count
    return isoform_count, splice_isoform_count, len(total_uniq_utrs), len(read_subsample)


# outfile = open("../../results/scratch/countCombinatorialIsoforms/combinatorial_isoform_count_subsampled.txt",'w')
# utr_assignments = {}
# total_utrs = open("../../results/utrs/assignments/all_isoforms_utrs.tsv")
#
# for line in total_utrs:
#     fields = line.strip().split()
#     read_id = fields[0]
#     # print line
#     gene_id = fields[1].split("-cluster")[0]
#     # cluster_id = int(fields[1].split("-cluster")[1])
#     cluster_id = fields[1]
#     # total_uniq_utrs.add(cluster_id)
#     utr_assignments[read_id] = (gene_id,cluster_id)
#
# outfile.write("stage\tdataset\tcounts\tfraction_retained\tread_count\n")
#
# total_in = "../../results/isoforms/all_isoforms.tsv"
# l1_in = "../../results/isoforms/L1_isoforms.tsv"
# l2_in = "../../results/isoforms/L2_isoforms.tsv"
# l3_in = "../../results/isoforms/L3_isoforms.tsv"
# l4_in = "../../results/isoforms/L4_isoforms.tsv"
# ya_in = "../../results/isoforms/YA_isoforms.tsv"
# ga_in = "../../results/isoforms/GA_isoforms.tsv"
# ml_in = "../../results/isoforms/ML_isoforms.tsv"
#
#
# for frac in np.linspace(0.001,1.0,10):
#     for x in range(1):
#     #for x in range(10):
#         print "%f - iter %d" %(frac,x)
#         all_isoform_count,all_splice_count,all_utr_count,all_read_count = countCombinatorialIsoforms(total_in,utr_assignments,fraction_retained=frac)
#         outfile.write("all\tfull length isoforms\t%d\t%f\t%d\n"%(all_isoform_count,frac,all_read_count))
#
#         l1_isoform_count,l1_splice_count,l1_utr_count,l1_read_count = countCombinatorialIsoforms(l1_in,utr_assignments,fraction_retained=frac)
#         outfile.write("L1\tfull length isoforms\t%d\t%f\t%d\n"%(l1_isoform_count,frac,l1_read_count))
#
#         l2_isoform_count,l2_splice_count,l2_utr_count,l2_read_count = countCombinatorialIsoforms(l2_in,utr_assignments,fraction_retained=frac)
#         outfile.write("L2\tfull length isoforms\t%d\t%f\t%d\n"%(l2_isoform_count,frac,l2_read_count))
#
#         l3_isoform_count,l3_splice_count,l3_utr_count,l3_read_count = countCombinatorialIsoforms(l3_in,utr_assignments,fraction_retained=frac)
#         outfile.write("L3\tfull length isoforms\t%d\t%f\t%d\n"%(l3_isoform_count,frac,l3_read_count))
#
#         l4_isoform_count,l4_splice_count,l4_utr_count,l4_read_count = countCombinatorialIsoforms(l4_in,utr_assignments,fraction_retained=frac)
#         outfile.write("L4\tfull length isoforms\t%d\t%f\t%d\n"%(l4_isoform_count,frac,l4_read_count))
#
#         ya_isoform_count,ya_splice_count,ya_utr_count,ya_read_count = countCombinatorialIsoforms(ya_in,utr_assignments,fraction_retained=frac)
#         outfile.write("young adult\tfull length isoforms\t%d\t%f\t%d\n"%(ya_isoform_count,frac,ya_read_count))
#
#         ga_isoform_count,ga_splice_count,ga_utr_count,ga_read_count = countCombinatorialIsoforms(ga_in,utr_assignments,fraction_retained=frac)
#         outfile.write("mature adult\tfull length isoforms\t%d\t%f\t%d\n"%(ga_isoform_count,frac,ga_read_count))
#
#         ml_isoform_count,ml_splice_count,ml_utr_count,ml_read_count = countCombinatorialIsoforms(ml_in,utr_assignments,fraction_retained=frac)
#         outfile.write("male\tfull length isoforms\t%d\t%f\t%d\n"%(ml_isoform_count,frac,ml_read_count))
#
#         outfile.flush()

# ### Sensitive
# outfile = open("../../results/scratch/countCombinatorialIsoforms/sensitive_combinatorial_isoform_count_subsampled.txt",'w')
# utr_assignments = {}
# total_utrs = open("../../results/utrs/assignments/all_sensitive_isoforms_utrs.tsv")
#
# for line in total_utrs:
#     fields = line.strip().split()
#     read_id = fields[0]
#     # print line
#     gene_id = fields[1].split("-cluster")[0]
#     # cluster_id = int(fields[1].split("-cluster")[1])
#     cluster_id = fields[1]
#     # total_uniq_utrs.add(cluster_id)
#     utr_assignments[read_id] = (gene_id,cluster_id)
#
# outfile.write("stage\tdataset\tcounts\tfraction_retained\tread_count\n")
#
# total_in = "../../results/isoforms/all_sensitive_isoforms.tsv"
# l1_in = "../../results/isoforms/L1_sensitive_isoforms.tsv"
# l2_in = "../../results/isoforms/L2_sensitive_isoforms.tsv"
# l3_in = "../../results/isoforms/L3_sensitive_isoforms.tsv"
# l4_in = "../../results/isoforms/L4_sensitive_isoforms.tsv"
# ya_in = "../../results/isoforms/YA_sensitive_isoforms.tsv"
# ga_in = "../../results/isoforms/GA_sensitive_isoforms.tsv"
# ml_in = "../../results/isoforms/ML_sensitive_isoforms.tsv"
#
#
# for frac in np.linspace(0.001,1.0,10):
#     for x in range(1):
#     #for x in range(10):
#         print "%f - iter %d" %(frac,x)
#         all_isoform_count,all_splice_count,all_utr_count,all_read_count = countCombinatorialIsoforms(total_in,utr_assignments,fraction_retained=frac)
#         outfile.write("all\tfull length isoforms\t%d\t%f\t%d\n"%(all_isoform_count,frac,all_read_count))
#
#         l1_isoform_count,l1_splice_count,l1_utr_count,l1_read_count = countCombinatorialIsoforms(l1_in,utr_assignments,fraction_retained=frac)
#         outfile.write("L1\tfull length isoforms\t%d\t%f\t%d\n"%(l1_isoform_count,frac,l1_read_count))
#
#         l2_isoform_count,l2_splice_count,l2_utr_count,l2_read_count = countCombinatorialIsoforms(l2_in,utr_assignments,fraction_retained=frac)
#         outfile.write("L2\tfull length isoforms\t%d\t%f\t%d\n"%(l2_isoform_count,frac,l2_read_count))
#
#         l3_isoform_count,l3_splice_count,l3_utr_count,l3_read_count = countCombinatorialIsoforms(l3_in,utr_assignments,fraction_retained=frac)
#         outfile.write("L3\tfull length isoforms\t%d\t%f\t%d\n"%(l3_isoform_count,frac,l3_read_count))
#
#         l4_isoform_count,l4_splice_count,l4_utr_count,l4_read_count = countCombinatorialIsoforms(l4_in,utr_assignments,fraction_retained=frac)
#         outfile.write("L4\tfull length isoforms\t%d\t%f\t%d\n"%(l4_isoform_count,frac,l4_read_count))
#
#         ya_isoform_count,ya_splice_count,ya_utr_count,ya_read_count = countCombinatorialIsoforms(ya_in,utr_assignments,fraction_retained=frac)
#         outfile.write("young adult\tfull length isoforms\t%d\t%f\t%d\n"%(ya_isoform_count,frac,ya_read_count))
#
#         ga_isoform_count,ga_splice_count,ga_utr_count,ga_read_count = countCombinatorialIsoforms(ga_in,utr_assignments,fraction_retained=frac)
#         outfile.write("mature adult\tfull length isoforms\t%d\t%f\t%d\n"%(ga_isoform_count,frac,ga_read_count))
#
#         ml_isoform_count,ml_splice_count,ml_utr_count,ml_read_count = countCombinatorialIsoforms(ml_in,utr_assignments,fraction_retained=frac)
#         outfile.write("male\tfull length isoforms\t%d\t%f\t%d\n"%(ml_isoform_count,frac,ml_read_count))
#
#         outfile.flush()



### Stringent
outfile = open("../../results/scratch/countCombinatorialIsoforms/stringent_combinatorial_isoform_count_subsampled.txt",'w')
utr_assignments = {}
total_utrs = open("../../results/utrs/assignments/all_stringent_isoforms_utrs.tsv")

for line in total_utrs:
    fields = line.strip().split()
    read_id = fields[0]
    # print line
    gene_id = fields[1].split("-cluster")[0]
    # cluster_id = int(fields[1].split("-cluster")[1])
    cluster_id = fields[1]
    # total_uniq_utrs.add(cluster_id)
    utr_assignments[read_id] = (gene_id,cluster_id)

outfile.write("stage\tdataset\tcounts\tfraction_retained\tread_count\n")

total_in = "../../results/isoforms/all_stringent_isoforms.tsv"
l1_in = "../../results/isoforms/L1_stringent_isoforms.tsv"
l2_in = "../../results/isoforms/L2_stringent_isoforms.tsv"
l3_in = "../../results/isoforms/L3_stringent_isoforms.tsv"
l4_in = "../../results/isoforms/L4_stringent_isoforms.tsv"
ya_in = "../../results/isoforms/YA_stringent_isoforms.tsv"
ga_in = "../../results/isoforms/GA_stringent_isoforms.tsv"
ml_in = "../../results/isoforms/ML_stringent_isoforms.tsv"


for frac in np.linspace(0.001,1.0,10):
    for x in range(1):
    #for x in range(10):
        print "%f - iter %d" %(frac,x)
        all_isoform_count,all_splice_count,all_utr_count,all_read_count = countCombinatorialIsoforms(total_in,utr_assignments,fraction_retained=frac)
        outfile.write("all\tfull length isoforms\t%d\t%f\t%d\n"%(all_isoform_count,frac,all_read_count))

        l1_isoform_count,l1_splice_count,l1_utr_count,l1_read_count = countCombinatorialIsoforms(l1_in,utr_assignments,fraction_retained=frac)
        outfile.write("L1\tfull length isoforms\t%d\t%f\t%d\n"%(l1_isoform_count,frac,l1_read_count))

        l2_isoform_count,l2_splice_count,l2_utr_count,l2_read_count = countCombinatorialIsoforms(l2_in,utr_assignments,fraction_retained=frac)
        outfile.write("L2\tfull length isoforms\t%d\t%f\t%d\n"%(l2_isoform_count,frac,l2_read_count))

        l3_isoform_count,l3_splice_count,l3_utr_count,l3_read_count = countCombinatorialIsoforms(l3_in,utr_assignments,fraction_retained=frac)
        outfile.write("L3\tfull length isoforms\t%d\t%f\t%d\n"%(l3_isoform_count,frac,l3_read_count))

        l4_isoform_count,l4_splice_count,l4_utr_count,l4_read_count = countCombinatorialIsoforms(l4_in,utr_assignments,fraction_retained=frac)
        outfile.write("L4\tfull length isoforms\t%d\t%f\t%d\n"%(l4_isoform_count,frac,l4_read_count))

        ya_isoform_count,ya_splice_count,ya_utr_count,ya_read_count = countCombinatorialIsoforms(ya_in,utr_assignments,fraction_retained=frac)
        outfile.write("young adult\tfull length isoforms\t%d\t%f\t%d\n"%(ya_isoform_count,frac,ya_read_count))

        ga_isoform_count,ga_splice_count,ga_utr_count,ga_read_count = countCombinatorialIsoforms(ga_in,utr_assignments,fraction_retained=frac)
        outfile.write("mature adult\tfull length isoforms\t%d\t%f\t%d\n"%(ga_isoform_count,frac,ga_read_count))

        ml_isoform_count,ml_splice_count,ml_utr_count,ml_read_count = countCombinatorialIsoforms(ml_in,utr_assignments,fraction_retained=frac)
        outfile.write("male\tfull length isoforms\t%d\t%f\t%d\n"%(ml_isoform_count,frac,ml_read_count))

        outfile.flush()
