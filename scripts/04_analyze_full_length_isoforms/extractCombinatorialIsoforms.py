#!/usr/bin/env python2

import sets

### nroach2@jhu.edu
### John has asked me to quantify the total number of isoforms we see in each stage and across all stages, including variation in 3'UTRs
### "For now don't worry about 3'UTR clustering"

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

# total_utrs = open("../results/realigned_utrs/assignments/all_isoforms_utrs.tsv")
# total_bed = open("../results/realigned_utrs/beds/all_isoforms_utrs.bed")
# utr_assignments = {}
# utr_end_points = {}
# for line in total_utrs:
#     fields = line.strip().split()
#     read_id = fields[0]
#     gene_id = fields[1].split("-cluster")[0]
#     cluster_id = fields[1]
#     utr_assignments[read_id] = (gene_id,cluster_id)
#
# for line in total_bed:
#     fields = line.strip().split()
#     strand = fields[5]
#     cluster_id = fields[3]
#     if strand == '+':
#         end = int(fields[2])
#     elif strand == '-':
#         end = int(fields[1])
#     utr_end_points[cluster_id] = end
#
# total_in = open("../results/realigned_isoforms/all_isoforms.tsv")
# outfile = open("combinatorial_isoforms.bed",'w')
# printCombinatorialIsoforms(total_in,utr_assignments,utr_end_points,outfile)
#

### Stringent
total_utrs = open("../../results/utrs/assignments/all_stringent_isoforms_utrs.tsv")
total_bed = open("../../results/utrs/beds/all_stringent_isoforms_utrs.bed")
utr_assignments = {}
utr_end_points = {}
for line in total_utrs:
    fields = line.strip().split()
    read_id = fields[0]
    gene_id = fields[1].split("-cluster")[0]
    cluster_id = fields[1]
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

total_in = open("../../results/isoforms/all_stringent_isoforms.tsv")
outfile = open("../../results/isoforms/stringent_combinatorial_isoforms.bed",'w')
printCombinatorialIsoforms(total_in,utr_assignments,utr_end_points,outfile)

# ###Sensitive
# total_utrs = open("../../results/utrs/assignments/all_sensitive_isoforms_utrs.tsv")
# total_bed = open("../../results/utrs/beds/all_sensitive_isoforms_utrs.bed")
# utr_assignments = {}
# utr_end_points = {}
# for line in total_utrs:
#     fields = line.strip().split()
#     read_id = fields[0]
#     gene_id = fields[1].split("-cluster")[0]
#     cluster_id = fields[1]
#     utr_assignments[read_id] = (gene_id,cluster_id)
#
# for line in total_bed:
#     fields = line.strip().split()
#     strand = fields[5]
#     cluster_id = fields[3]
#     if strand == '+':
#         end = int(fields[2])
#     elif strand == '-':
#         end = int(fields[1])
#     utr_end_points[cluster_id] = end
#
# total_in = open("../../results/isoforms/all_sensitive_isoforms.tsv")
# outfile = open("../../results/isoforms/sensitive_combinatorial_isoforms.bed",'w')
# printCombinatorialIsoforms(total_in,utr_assignments,utr_end_points,outfile)