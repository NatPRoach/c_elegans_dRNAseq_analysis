#!/usr/bin/env python2

import sets

### nroach2@jhu.edu
### John has asked me to quantify the total number of isoforms we see in each stage and across all stages, including variation in 3'UTRs
### "For now don't worry about 3'UTR clustering"

# def getStartCodon(strand,start,block_sizes,block_starts,cds_start):
#     if strand == '+':
#         for x in range(len(block_starts)):
#             if cds_start > block_sizes[x]:
#                 cds_start -= block_sizes[x]
#             else:
#                 return start + block_starts[x] + cds_start
#     elif strand == '-':
#         for x in range(len(block_starts)):
#             if cds_start > block_sizes[len(block_starts)-1-x]:
#                 cds_start -= block_sizes[len(block_starts)-1-x]
#             else:
#                 return start + block_starts[len(block_starts)-1-x] + (block_sizes[len(block_starts)-1-x] - cds_start)
#
#
# def getStartPosition(cds_start,start_codon,introns,strand):
#     if strand == '+':
#
#     elif strand == '-':
#
#
#     return start

def convertIntronsToBlockStartsAndSizes(introns,start,end):
    sorted_introns = []
    for intron in introns:
        if intron[0] > intron[1]:
            new_intron = (intron[1],intron[0])
            sorted_introns.append(new_intron)
        else:
            sorted_introns.append(intron)
    sorted_introns.sort()
    
    block_sizes = []
    block_starts = []
    
    next_start = start
    for i1,i2 in sorted_introns:
        block_sizes.append(i1 - next_start)
        block_starts.append(next_start - start)
        next_start = i2
    block_sizes.append(end-next_start)
    block_starts.append(next_start - start)
    
    return block_sizes,block_starts

def printCombinatorialIsoforms(infile,utr_assignments,utr_end_points,outfile):
    isoform_count = 0
    splice_isoform_count = 0
    total_uniq_utrs = sets.Set()
    gene_spliceform_counter = {}
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
        if len(fields ) > 15:
            introns = [(int(x.split(',')[0]),int(x.split(',')[1])) for x in fields[15].split(';')]
        else:
            introns = []
        assert len(stop_positions) == len(read_ids)
        if retained_intron: #For now ignore retained intron transcripts
            continue
        if gene not in gene_spliceform_counter:
            gene_spliceform_counter[gene] = 0
        else:
            gene_spliceform_counter[gene] += 1
        #Each line represents an intron chain (with 5' tolerance... remove 5' tolerance?)
        #Need to count number of utrs represented in each intron chain
        #We have a file of read_id -> gene & UTR call
        uniq_utrs = sets.Set()
        utr_support_count = {}
        for read_id in read_ids:
            if read_id in utr_assignments:
                gene_id,cluster_id = utr_assignments[read_id]
                assert gene_id == gene
                uniq_utrs.add(cluster_id)
                if cluster_id in utr_support_count:
                    utr_support_count[cluster_id] += 1
                else:
                    utr_support_count[cluster_id] = 1
                total_uniq_utrs.add(cluster_id)
        if strand == '+':
            start_point = min(start_positions)
        elif strand == '-':
            start_point = max(start_positions)
        if len(uniq_utrs) == 0:
            #No defined utrs here, use the longest utr
            isoform_label = gene + "-Spliceform" + str(gene_spliceform_counter[gene]) + "-NoUTR"
            if strand == '+':
                end_point = max(stop_positions)
                block_sizes,block_starts = convertIntronsToBlockStartsAndSizes(introns,start_point,end_point)
                outfile.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t255,0,0\t%d\t%s\t%s\n" %(chrom,start_point,end_point,isoform_label,min((1000,num_reads)),strand,start_codon,stop_codon,len(block_sizes),','.join([str(x) for x in block_sizes]),','.join([str(x) for x in block_starts])))
            elif strand == '-':
                end_point = min(stop_positions)
                block_sizes,block_starts = convertIntronsToBlockStartsAndSizes(introns,end_point,start_point)
                outfile.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t255,0,0\t%d\t%s\t%s\n" %(chrom,end_point,start_point,isoform_label,min((1000,num_reads)),strand,stop_codon,start_codon,len(block_sizes),','.join([str(x) for x in block_sizes]),','.join([str(x) for x in block_starts])))
            
        else:
            for cluster_id in uniq_utrs:
                utr_id = cluster_id.split("-cluster")[1]
                isoform_label = gene + "-Spliceform" + str(gene_spliceform_counter[gene]) + "-UTR" + utr_id
                end_point = utr_end_points[cluster_id]
                if strand == '+':
                    block_sizes,block_starts = convertIntronsToBlockStartsAndSizes(introns,start_point,end_point)
                    outfile.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t0,0,255\t%d\t%s\t%s\n" %(chrom,start_point,end_point,isoform_label,min((1000,utr_support_count[cluster_id])),strand,start_codon,stop_codon,len(block_sizes),','.join([str(x) for x in block_sizes]),','.join([str(x) for x in block_starts])))
                elif strand == '-':
                    block_sizes,block_starts = convertIntronsToBlockStartsAndSizes(introns,end_point,start_point)
                    outfile.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t0,0,255\t%d\t%s\t%s\n" %(chrom,end_point,start_point,isoform_label,min((1000,utr_support_count[cluster_id])),strand,stop_codon,start_codon,len(block_sizes),','.join([str(x) for x in block_sizes]),','.join([str(x) for x in block_starts])))

# def countCombinatorialIsoforms(infile,utr_assignments,utr_end_points,outfile):
#     isoform_count = 0
#     splice_isoform_count = 0
#     total_uniq_utrs = sets.Set()
#     for line in infile:
#         splice_isoform_count += 1
#         fields = line.strip().split()
#         gene = fields[0]
#         chrom = fields[1]
#         strand = fields[2]
#         cds = fields[3]
#         seq = fields[4]
#         cds_start = int(fields[5])
#         cds_end = int(fields[6])
#         num_reads = int(fields[7])
#         num_exons = int(fields[8])
#         start_codon = int(fields[9])
#         stop_codon = int(fields[10])
#         stop_positions = [int(x) for x in fields[11].split(',')]
#         read_ids = fields[12].split(',')
#         retained_intron = fields[13] == '1'
#         assert len(stop_positions) == len(read_ids)
#         if retained_intron: #For now ignore retained intron transcripts
#             continue
#         #Each line represents an intron chain (with 5' tolerance... remove 5' tolerance?)
#         #Need to count number of utrs represented in each intron chain
#         #We have a file of read_id -> gene & UTR call
#         uniq_utrs = sets.Set()
#         for read_id in read_ids:
#             if read_id in utr_assignments:
#                 gene_id,cluster_id = utr_assignments[read_id]
#                 assert gene_id == gene
#                 uniq_utrs.add(cluster_id)
#                 total_uniq_utrs.add(cluster_id)
#         if len(uniq_utrs) == 0:
#             utr_count = 1
#         else:
#             utr_count = len(uniq_utrs)
#         isoform_count += utr_count
#     return isoform_count, splice_isoform_count, len(total_uniq_utrs)


total_utrs = open("../results/realigned_utrs/assignments/all_isoforms_utrs.tsv")
total_bed = open("../results/realigned_utrs/beds/all_isoforms_utrs.bed")
utr_assignments = {}
# total_uniq_utrs = sets.Set()
utr_end_points = {}
for line in total_utrs:
    fields = line.strip().split()
    read_id = fields[0]
    gene_id = fields[1].split("-cluster")[0]
    # cluster_id = int(fields[1].split("-cluster")[1])
    cluster_id = fields[1]
    # total_uniq_utrs.add(cluster_id)
    utr_assignments[read_id] = (gene_id,cluster_id)

for line in total_bed:
    fields = line.strip().split()
    strand = fields[5]
    cluster_id = fields[3]
    if strand == '+':
        end = int(fields[2])
    elif strand == '-':
        end = int(fields[1])
    utr_end_points[cluster_id] = end

# for uniq_utr in total_uniq_utrs:
#     print uniq_utr
# print total_uniq_utrs
# print len(total_uniq_utrs)

total_in = open("../results/realigned_isoforms/all_isoforms.tsv")
outfile = open("combinatorial_isoforms.bed",'w')
printCombinatorialIsoforms(total_in,utr_assignments,utr_end_points,outfile)
# l1_in = open("../results/realigned_isoforms/L1_isoforms.tsv")
# outfile = open("combinatorial_isoforms.bed",'w')
# printCombinatorialIsoforms(l1_in,utr_assignments,utr_end_points,outfile)


#all_isoform_count,all_splice_count,all_utr_count = countCombinatorialIsoforms(total_in,utr_assignments,utr_end_points,outfile)
#
# l1_in = open("../results/realigned_isoforms/L1_isoforms.tsv")
# l1_isoform_count,l1_splice_count,l1_utr_count = countCombinatorialIsoforms(l1_in,utr_assignments)
#
# l2_in = open("../results/realigned_isoforms/L2_isoforms.tsv")
# l2_isoform_count,l2_splice_count,l2_utr_count = countCombinatorialIsoforms(l2_in,utr_assignments)
#
# l3_in = open("../results/realigned_isoforms/L3_isoforms.tsv")
# l3_isoform_count,l3_splice_count,l3_utr_count = countCombinatorialIsoforms(l3_in,utr_assignments)
#
# l4_in = open("../results/realigned_isoforms/L4_isoforms.tsv")
# l4_isoform_count,l4_splice_count,l4_utr_count = countCombinatorialIsoforms(l4_in,utr_assignments)
#
# ya_in = open("../results/realigned_isoforms/YA_isoforms.tsv")
# ya_isoform_count,ya_splice_count,ya_utr_count = countCombinatorialIsoforms(ya_in,utr_assignments)
#
# ga_in = open("../results/realigned_isoforms/GA_isoforms.tsv")
# ga_isoform_count,ga_splice_count,ga_utr_count = countCombinatorialIsoforms(ga_in,utr_assignments)
#
# ml_in = open("../results/realigned_isoforms/ML_isoforms.tsv")
# ml_isoform_count,ml_splice_count,ml_utr_count = countCombinatorialIsoforms(ml_in,utr_assignments)
