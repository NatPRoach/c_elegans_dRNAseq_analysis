#!/usr/bin/env python2

from __future__ import division

from matplotlib import pyplot as plt
import numpy as np
# import pymc3 as pm
# import scipy as sp
import seaborn as sns
# from theano import tensor as tt
# import pandas as pd
import time
import sys
import sets

def calculatePutativeOffsets(values,scaling_factor=0.0,kernal_width = 10,window = 10,proportion_threshold=0.0005,read_threshold=3):
    bw = max(np.std(values)*scaling_factor,kernal_width)
    ax = sns.kdeplot(values,bw = bw)
    x_vals,y_vals = ax.lines[0].get_data()
    counter = 0
    x_pos = []
    for y in range(1,len(y_vals)-1):
        x_counter = 0
        for val in values:
            if abs(x_vals[y] - val) < window:
                x_counter += 1
        if x_counter < read_threshold:
            continue
        if y_vals[y-1] < y_vals[y] and y_vals[y+1] < y_vals[y]:
            x_pos.append(x_vals[y])
    plt.clf()
    return x_pos
def model(values,strand = '+',read_threshold=3):
    new_values = []
    if strand == '+':
        max_value = max(values)
        for x in values:
            new_values.append(-1*(x - max_value))
    elif strand == '-':
        min_value = min(values)
        for x in values:
            new_values.append(x-min_value)
    
    putative_offsets = calculatePutativeOffsets(new_values,read_threshold=3)
    endpoints2 = []
    if strand == '+':
        for endpoint in putative_offsets:
            endpoints2.append(int(round(max_value - endpoint)))
    elif strand == '-':
        for endpoint in putative_offsets:
            endpoints2.append(int(round(min_value + endpoint)))
    # return [endpoints, endpoints2]
    return endpoints2

def convertIntronsToBedBlocks(stop_codon,end_point,introns_past_stop,strand):
    block_starts = []
    block_sizes  = []
    if strand == '+':
        #total_len = end_point - stop_codon
        last_block_start = stop_codon 
        for intron in introns_past_stop:
            block_size = intron[0] - stop_codon
            block_start = last_block_start - stop_codon
            block_sizes.append(block_size)
            block_starts.append(block_start)
            last_block_start = intron[1]
        block_sizes.append(end_point - last_block_start)
        block_starts.append(last_block_start - stop_codon)
    elif strand == '-':
        #total_len = stop_codon - end_point
        last_block_start = end_point
        for intron in reversed(introns_past_stop):
            intron2 = (intron[1],intron[0])
            block_size = intron2[0] - end_point
            block_start = last_block_start - end_point
            block_sizes.append(block_size)
            block_starts.append(block_start)
            last_block_start = intron2[1]
        block_sizes.append(stop_codon - last_block_start + 1 )
        block_starts.append(last_block_start - end_point)
    return [block_starts,block_sizes]

plt.rcParams['figure.figsize'] = 14, 6
SEED = 4662 #random.org 
np.random.seed(SEED)
## Set basic parameters:
window = 10
minEndPoints = 3
maxUTRLength = 10000
percentage = 0.10
read_threshold=3
last_gene = None
last_chrom = None
last_strand = None
stop_position_dict = {}
read_id_dict = {}
infile = open(sys.argv[1])
outfile = open(sys.argv[2],'w')
assignment_file = open(sys.argv[3],'w')
exclusion_reads_in = open("../../data/scratch/all.exclude.txt")

all_clusters = sets.Set()
all_clusters2= sets.Set()
flag = False
cluster_sizes = []
counter = 0
# total_counter = 0
exclusion_reads = sets.Set()
for line in exclusion_reads_in:
    read_id = line.strip()
    exclusion_reads.add(read_id)
    
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
    start_positions = [int(x) for x in fields[11].split(',')]
    stop_positions = [int(x) for x in fields[12].split(',')]
    read_ids = fields[13].split(',')
    retained_intron = fields[14] == '1'
    if len(fields) == 16:
        introns = []
        for x in fields[15].split(';'):
            y,z = x.split(',')
            introns.append((int(y),int(z)))
    
        introns_past_stop = []
        for intron in introns:
            if strand == '+':
                if intron[0] > stop_codon:
                    introns_past_stop.append(intron)
            elif strand == '-':
                if intron[0] < stop_codon:
                    introns_past_stop.append(intron)
        introns_past_stop = tuple(introns_past_stop)
    else:
        introns_past_stop = ()
    assert len(stop_positions) == len(read_ids)
    if retained_intron: #For now ignore retained intron transcripts
        continue
    if len(read_ids) == 0:
        continue
    if gene == last_gene:
        if (stop_codon, introns_past_stop) not in stop_position_dict:
            stop_position_dict[(stop_codon,introns_past_stop)] = []
            read_id_dict[(stop_codon,introns_past_stop)] = []
        for x in range(len(stop_positions)):
            read_id = read_ids[x]
            stop_pos = stop_positions[x]
            if read_id not in exclusion_reads:
                if (strand == '+' and stop_pos > stop_codon) or (strand == '-' and stop_pos < stop_codon):
                    stop_position_dict[(stop_codon,introns_past_stop)].append(stop_pos)
                    read_id_dict[(stop_codon,introns_past_stop)].append(read_id)
    else:
        for stop_codon2,introns_past_stop2 in stop_position_dict:
            previous_cluster_count = counter 
            stop_pos = stop_position_dict[(stop_codon2,introns_past_stop2)]
            read_ids2 = read_id_dict[(stop_codon2,introns_past_stop2)]
            if len(stop_pos) >= read_threshold:
                end_points = model(stop_pos,strand = last_strand,read_threshold = read_threshold)
                cluster_sizes.append(len(end_points))
                for end_point in end_points:
                    if introns_past_stop2 != ():
                        block_starts, block_sizes = convertIntronsToBedBlocks(stop_codon2,end_point,introns_past_stop2,last_strand)
                    else:
                        block_starts = [0]
                        if last_strand == '-':
                            block_sizes = [abs(stop_codon2 - int(end_point))]
                        elif last_strand == '+':
                            block_sizes = [abs(stop_codon2 - int(end_point))]
                    if last_strand == '+':
                        outfile.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t0,0,255\t%d\t%s\t%s\n" %(last_chrom,stop_codon2,end_point,last_gene + "-cluster%d" %(counter),1000,last_strand,stop_codon2,stop_codon2,len(block_starts),','.join([str(x) for x in block_sizes]),','.join([str(x) for x in block_starts])))
                    elif last_strand == '-':
                        outfile.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t0,0,255\t%d\t%s\t%s\n" %(last_chrom,end_point,stop_codon2,last_gene + "-cluster%d" %(counter),1000,last_strand,end_point,end_point,len(block_starts),','.join([str(x) for x in block_sizes]),','.join([str(x) for x in block_starts])))
                    counter += 1
                for x in range(len(stop_pos)):
                    read_id = read_ids2[x]
                    stop = stop_pos[x]
                    
                    min_dist = np.inf
                    min_idx = None
                    for i,y in enumerate(end_points):
                        flag = True
                        if abs(y-stop) < min_dist:
                            if abs(y - stop) < 10.:
                                min_dist = abs(y - stop)
                                min_idx = i
                    if min_idx is not None:
                        assignment_file.write("%s\t%s\n" %(read_id,last_gene + "-cluster%d" %(min_idx + previous_cluster_count)))
        last_gene = gene
        last_chrom = chrom
        last_strand = strand
        stop_position_dict = {}
        read_id_dict = {}
        counter = 0
        stop_position_dict[(stop_codon,introns_past_stop)] = []
        read_id_dict[(stop_codon,introns_past_stop)] = []
        for x in range(len(stop_positions)):
            read_id = read_ids[x]
            stop_pos = stop_positions[x]
            if read_id not in exclusion_reads:
                if (strand == '+' and stop_pos > stop_codon) or (strand == '-' and stop_pos < stop_codon):
                    stop_position_dict[(stop_codon,introns_past_stop)].append(stop_pos)
                    read_id_dict[(stop_codon,introns_past_stop)].append(read_id)