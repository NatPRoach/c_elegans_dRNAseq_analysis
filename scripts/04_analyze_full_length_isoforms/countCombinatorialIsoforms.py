#!/usr/bin/env python2

import sets

### nroach2@jhu.edu
### John has asked me to quantify the total number of isoforms we see in each stage and across all stages, including variation in 3'UTRs
### "For now don't worry about 3'UTR clustering"

def countCombinatorialIsoforms(infile,utr_assignments):
    isoform_count = 0
    splice_isoform_count = 0
    total_uniq_utrs = sets.Set()
    for line in infile:
        splice_isoform_count += 1
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
        start_positions = [int(x) for x in fields[11].split(',')]
        stop_positions = [int(x) for x in fields[12].split(',')]
        read_ids = fields[13].split(',')
        retained_intron = fields[14] == '1'
        assert len(stop_positions) == len(read_ids)
        if retained_intron: #For now ignore retained intron transcripts
            continue
        #Each line represents an intron chain (with 5' tolerance... remove 5' tolerance?)
        #Need to count number of utrs represented in each intron chain
        #We have a file of read_id -> gene & UTR call
        uniq_utrs = sets.Set()
        for read_id in read_ids:
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
    return isoform_count, splice_isoform_count, len(total_uniq_utrs)


# total_utrs = open("../../results/utrs/assignments/all_isoforms_utrs.tsv")
# outfile = open("../../results/scratch/countCombinatorialIsoforms/combinatorial_isoform_count.txt",'w')
# utr_assignments = {}
# for line in total_utrs:
#     fields = line.strip().split()
#     read_id = fields[0]
#     gene_id = fields[1].split("-cluster")[0]\
#     cluster_id = fields[1]
#     utr_assignments[read_id] = (gene_id,cluster_id)
#
# total_in = open("../../results/isoforms/all_isoforms.tsv")
# all_isoform_count,all_splice_count,all_utr_count = countCombinatorialIsoforms(total_in,utr_assignments)
#
# l1_in = open("../../results/isoforms/L1_isoforms.tsv")
# l1_isoform_count,l1_splice_count,l1_utr_count = countCombinatorialIsoforms(l1_in,utr_assignments)
#
# l2_in = open("../../results/isoforms/L2_isoforms.tsv")
# l2_isoform_count,l2_splice_count,l2_utr_count = countCombinatorialIsoforms(l2_in,utr_assignments)
#
# l3_in = open("../../results/isoforms/L3_isoforms.tsv")
# l3_isoform_count,l3_splice_count,l3_utr_count = countCombinatorialIsoforms(l3_in,utr_assignments)
#
# l4_in = open("../../results/isoforms/L4_isoforms.tsv")
# l4_isoform_count,l4_splice_count,l4_utr_count = countCombinatorialIsoforms(l4_in,utr_assignments)
#
# ya_in = open("../../results/isoforms/YA_isoforms.tsv")
# ya_isoform_count,ya_splice_count,ya_utr_count = countCombinatorialIsoforms(ya_in,utr_assignments)
#
# ga_in = open("../../results/isoforms/GA_isoforms.tsv")
# ga_isoform_count,ga_splice_count,ga_utr_count = countCombinatorialIsoforms(ga_in,utr_assignments)
#
# ml_in = open("../../results/isoforms/ML_isoforms.tsv")
# ml_isoform_count,ml_splice_count,ml_utr_count = countCombinatorialIsoforms(ml_in,utr_assignments)
#
# outfile.write("stage\tdataset\tcounts\tx_order\n")
#
# outfile.write("L1\tsplice isoforms\t%d\t0\n"%(l1_splice_count))
# outfile.write("L2\tsplice isoforms\t%d\t1\n"%(l2_splice_count))
# outfile.write("L3\tsplice isoforms\t%d\t2\n"%(l3_splice_count))
# outfile.write("L4\tsplice isoforms\t%d\t3\n"%(l4_splice_count))
# outfile.write("young adult\tsplice isoforms\t%d\t4\n"%(ya_splice_count))
# outfile.write("mature adult\tsplice isoforms\t%d\t5\n"%(ga_splice_count))
# outfile.write("male\tsplice isoforms\t%d\t6\n"%(ml_splice_count))
# outfile.write("all\tsplice isoforms\t%d\t7\n"%(all_splice_count))
#
# outfile.write("L1\tutrs\t%d\t0\n"%(l1_utr_count))
# outfile.write("L2\tutrs\t%d\t1\n"%(l2_utr_count))
# outfile.write("L3\tutrs\t%d\t2\n"%(l3_utr_count))
# outfile.write("L4\tutrs\t%d\t3\n"%(l4_utr_count))
# outfile.write("young adult\tutrs\t%d\t4\n"%(ya_utr_count))
# outfile.write("mature adult\tutrs\t%d\t5\n"%(ga_utr_count))
# outfile.write("male\tutrs\t%d\t6\n"%(ml_utr_count))
# outfile.write("all\tutrs\t%d\t7\n"%(all_utr_count))
#
# outfile.write("L1\tfull length isoforms\t%d\t0\n"%(l1_isoform_count))
# outfile.write("L2\tfull length isoforms\t%d\t1\n"%(l2_isoform_count))
# outfile.write("L3\tfull length isoforms\t%d\t2\n"%(l3_isoform_count))
# outfile.write("L4\tfull length isoforms\t%d\t3\n"%(l4_isoform_count))
# outfile.write("young adult\tfull length isoforms\t%d\t4\n"%(ya_isoform_count))
# outfile.write("mature adult\tfull length isoforms\t%d\t5\n"%(ga_isoform_count))
# outfile.write("male\tfull length isoforms\t%d\t6\n"%(ml_isoform_count))
# outfile.write("all\tfull length isoforms\t%d\t7\n"%(all_isoform_count))

# ### Sensitive
# total_utrs = open("../../results/utrs/assignments/all_sensitive_isoforms_utrs.tsv")
# outfile = open("../../results/scratch/countCombinatorialIsoforms/sensitive_combinatorial_isoform_count.txt",'w')
# utr_assignments = {}
# for line in total_utrs:
#     fields = line.strip().split()
#     read_id = fields[0]
#     gene_id = fields[1].split("-cluster")[0]
#     cluster_id = fields[1]
#     utr_assignments[read_id] = (gene_id,cluster_id)
#
# total_in = open("../../results/isoforms/all_sensitive_isoforms.tsv")
# all_isoform_count,all_splice_count,all_utr_count = countCombinatorialIsoforms(total_in,utr_assignments)
#
# l1_in = open("../../results/isoforms/L1_sensitive_isoforms.tsv")
# l1_isoform_count,l1_splice_count,l1_utr_count = countCombinatorialIsoforms(l1_in,utr_assignments)
#
# l2_in = open("../../results/isoforms/L2_sensitive_isoforms.tsv")
# l2_isoform_count,l2_splice_count,l2_utr_count = countCombinatorialIsoforms(l2_in,utr_assignments)
#
# l3_in = open("../../results/isoforms/L3_sensitive_isoforms.tsv")
# l3_isoform_count,l3_splice_count,l3_utr_count = countCombinatorialIsoforms(l3_in,utr_assignments)
#
# l4_in = open("../../results/isoforms/L4_sensitive_isoforms.tsv")
# l4_isoform_count,l4_splice_count,l4_utr_count = countCombinatorialIsoforms(l4_in,utr_assignments)
#
# ya_in = open("../../results/isoforms/YA_sensitive_isoforms.tsv")
# ya_isoform_count,ya_splice_count,ya_utr_count = countCombinatorialIsoforms(ya_in,utr_assignments)
#
# ga_in = open("../../results/isoforms/GA_sensitive_isoforms.tsv")
# ga_isoform_count,ga_splice_count,ga_utr_count = countCombinatorialIsoforms(ga_in,utr_assignments)
#
# ml_in = open("../../results/isoforms/ML_sensitive_isoforms.tsv")
# ml_isoform_count,ml_splice_count,ml_utr_count = countCombinatorialIsoforms(ml_in,utr_assignments)
#
# outfile.write("stage\tdataset\tcounts\tx_order\n")
#
# outfile.write("L1\tsplice isoforms\t%d\t0\n"%(l1_splice_count))
# outfile.write("L2\tsplice isoforms\t%d\t1\n"%(l2_splice_count))
# outfile.write("L3\tsplice isoforms\t%d\t2\n"%(l3_splice_count))
# outfile.write("L4\tsplice isoforms\t%d\t3\n"%(l4_splice_count))
# outfile.write("young adult\tsplice isoforms\t%d\t4\n"%(ya_splice_count))
# outfile.write("mature adult\tsplice isoforms\t%d\t5\n"%(ga_splice_count))
# outfile.write("male\tsplice isoforms\t%d\t6\n"%(ml_splice_count))
# outfile.write("all\tsplice isoforms\t%d\t7\n"%(all_splice_count))
#
# outfile.write("L1\tutrs\t%d\t0\n"%(l1_utr_count))
# outfile.write("L2\tutrs\t%d\t1\n"%(l2_utr_count))
# outfile.write("L3\tutrs\t%d\t2\n"%(l3_utr_count))
# outfile.write("L4\tutrs\t%d\t3\n"%(l4_utr_count))
# outfile.write("young adult\tutrs\t%d\t4\n"%(ya_utr_count))
# outfile.write("mature adult\tutrs\t%d\t5\n"%(ga_utr_count))
# outfile.write("male\tutrs\t%d\t6\n"%(ml_utr_count))
# outfile.write("all\tutrs\t%d\t7\n"%(all_utr_count))
#
# outfile.write("L1\tfull length isoforms\t%d\t0\n"%(l1_isoform_count))
# outfile.write("L2\tfull length isoforms\t%d\t1\n"%(l2_isoform_count))
# outfile.write("L3\tfull length isoforms\t%d\t2\n"%(l3_isoform_count))
# outfile.write("L4\tfull length isoforms\t%d\t3\n"%(l4_isoform_count))
# outfile.write("young adult\tfull length isoforms\t%d\t4\n"%(ya_isoform_count))
# outfile.write("mature adult\tfull length isoforms\t%d\t5\n"%(ga_isoform_count))
# outfile.write("male\tfull length isoforms\t%d\t6\n"%(ml_isoform_count))
# outfile.write("all\tfull length isoforms\t%d\t7\n"%(all_isoform_count))


### Stringent
total_utrs = open("../../results/utrs/assignments/all_stringent_isoforms_utrs.tsv")
outfile = open("../../results/scratch/countCombinatorialIsoforms/stringent_combinatorial_isoform_count.txt",'w')
utr_assignments = {}
for line in total_utrs:
    fields = line.strip().split()
    read_id = fields[0]
    gene_id = fields[1].split("-cluster")[0]
    cluster_id = fields[1]
    utr_assignments[read_id] = (gene_id,cluster_id)

total_in = open("../../results/isoforms/all_stringent_isoforms.tsv")
all_isoform_count,all_splice_count,all_utr_count = countCombinatorialIsoforms(total_in,utr_assignments)

l1_in = open("../../results/isoforms/L1_stringent_isoforms.tsv")
l1_isoform_count,l1_splice_count,l1_utr_count = countCombinatorialIsoforms(l1_in,utr_assignments)

l2_in = open("../../results/isoforms/L2_stringent_isoforms.tsv")
l2_isoform_count,l2_splice_count,l2_utr_count = countCombinatorialIsoforms(l2_in,utr_assignments)

l3_in = open("../../results/isoforms/L3_stringent_isoforms.tsv")
l3_isoform_count,l3_splice_count,l3_utr_count = countCombinatorialIsoforms(l3_in,utr_assignments)

l4_in = open("../../results/isoforms/L4_stringent_isoforms.tsv")
l4_isoform_count,l4_splice_count,l4_utr_count = countCombinatorialIsoforms(l4_in,utr_assignments)

ya_in = open("../../results/isoforms/YA_stringent_isoforms.tsv")
ya_isoform_count,ya_splice_count,ya_utr_count = countCombinatorialIsoforms(ya_in,utr_assignments)

ga_in = open("../../results/isoforms/GA_stringent_isoforms.tsv")
ga_isoform_count,ga_splice_count,ga_utr_count = countCombinatorialIsoforms(ga_in,utr_assignments)

ml_in = open("../../results/isoforms/ML_stringent_isoforms.tsv")
ml_isoform_count,ml_splice_count,ml_utr_count = countCombinatorialIsoforms(ml_in,utr_assignments)

outfile.write("stage\tdataset\tcounts\tx_order\n")

outfile.write("L1\tsplice isoforms\t%d\t0\n"%(l1_splice_count))
outfile.write("L2\tsplice isoforms\t%d\t1\n"%(l2_splice_count))
outfile.write("L3\tsplice isoforms\t%d\t2\n"%(l3_splice_count))
outfile.write("L4\tsplice isoforms\t%d\t3\n"%(l4_splice_count))
outfile.write("young adult\tsplice isoforms\t%d\t4\n"%(ya_splice_count))
outfile.write("mature adult\tsplice isoforms\t%d\t5\n"%(ga_splice_count))
outfile.write("male\tsplice isoforms\t%d\t6\n"%(ml_splice_count))
outfile.write("all\tsplice isoforms\t%d\t7\n"%(all_splice_count))

outfile.write("L1\tutrs\t%d\t0\n"%(l1_utr_count))
outfile.write("L2\tutrs\t%d\t1\n"%(l2_utr_count))
outfile.write("L3\tutrs\t%d\t2\n"%(l3_utr_count))
outfile.write("L4\tutrs\t%d\t3\n"%(l4_utr_count))
outfile.write("young adult\tutrs\t%d\t4\n"%(ya_utr_count))
outfile.write("mature adult\tutrs\t%d\t5\n"%(ga_utr_count))
outfile.write("male\tutrs\t%d\t6\n"%(ml_utr_count))
outfile.write("all\tutrs\t%d\t7\n"%(all_utr_count))

outfile.write("L1\tfull length isoforms\t%d\t0\n"%(l1_isoform_count))
outfile.write("L2\tfull length isoforms\t%d\t1\n"%(l2_isoform_count))
outfile.write("L3\tfull length isoforms\t%d\t2\n"%(l3_isoform_count))
outfile.write("L4\tfull length isoforms\t%d\t3\n"%(l4_isoform_count))
outfile.write("young adult\tfull length isoforms\t%d\t4\n"%(ya_isoform_count))
outfile.write("mature adult\tfull length isoforms\t%d\t5\n"%(ga_isoform_count))
outfile.write("male\tfull length isoforms\t%d\t6\n"%(ml_isoform_count))
outfile.write("all\tfull length isoforms\t%d\t7\n"%(all_isoform_count))
