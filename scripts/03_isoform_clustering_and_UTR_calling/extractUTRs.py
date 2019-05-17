#!/usr/bin/env python2

from __future__ import division

from matplotlib import pyplot as plt
import numpy as np
import pymc3 as pm
import scipy as sp
import seaborn as sns
from theano import tensor as tt
import pandas as pd
import time
import sys
import sets

def calculatePutativeOffsets(values,scaling_factor=0.0,kernal_width = 10,window = 10,proportion_threshold=0.0005,read_threshold=3):
    bw = max(np.std(values)*scaling_factor,kernal_width)
    #bw = kernal_width
    ax = sns.kdeplot(values,bw = bw)
    x_vals,y_vals = ax.lines[0].get_data()
    counter = 0
    x_pos = []
    for y in range(1,len(y_vals)-1):
        # if y_vals[y] <= proportion_threshold:
        #    continue
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
    # cluster_count = len(putative_offsets)
    # if cluster_count == 0:
    #     raw_endpoints = []
    #     return [[],[]]
    # elif cluster_count == 1:
    #     #return [[0],[0]]
    #     raw_endpoints = model1(new_values)
    # elif cluster_count == 2:
    #     #return [0,0]
    #     raw_endpoints = model2(new_values)
    # elif cluster_count == 3:
    #     raw_endpoints = model3(new_values)
    # elif cluster_count == 4:
    #     raw_endpoints = model4(new_values)
    # elif cluster_count == 5:
    #     raw_endpoints = model5(new_values)
    # elif cluster_count == 6:
    #     raw_endpoints = model6(new_values)
    # elif cluster_count == 7:
    #     raw_endpoints = model7(new_values)
    # elif cluster_count == 8:
    #     raw_endpoints = model8(new_values)
    # else:
    #     print ">8 clusters!, consider increasing kernel width or read threshold"
    #     raw_endpoints = []
    # endpoints = []
    # if strand == '+':
    #     for endpoint in raw_endpoints:
    #         endpoints.append(max_value - endpoint)
    # elif strand == '-':
    #     for endpoint in raw_endpoints:
    #         endpoints.append(min_value + endpoint)
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

#
# def model1(values):
#     return [0]
#     with pm.Model() as model:
#         # define priors, weakly informative Normal
#         b0 = pm.Normal('b0_intercept', mu=0, sd=10)
#         b1 = pm.Normal('b1_alcohol[T.True]', mu=0, sd=10)
#         #theta = (b0 + b1*)
#         # lam1 = pm.Exponential('lam1', lam=1.)
# #         offset1 = pm.Uniform('offset1',upper=float(max(values)))
#
#         pois1 = pm.Poisson.dist(mu=lam1+offset1,observed=values)
#         # like = pm.Model('like',w=w,comp_dists=[pois1],observed=values)
#     with model:
#         trace = pm.sample(50000,n_init=10000,tune=10000)
#     of1 = np.mean(trace["offset1"][-4000:])
#     lam1_avg = np.mean(trace['lam1'][-4000:])
#
#     cluster1Endpoint = int(round(of1 + lam1_avg))
#
#     return [cluster1Endpoint]
#
# def model2(values):
#     return [0 for x in range(2)]
#     with pm.Model() as model:
#         lam1 = pm.Exponential('lam1', lam=1.)
#         offset1 = pm.Uniform('offset1',upper=float(max(values)))
#         lam2 = pm.Exponential('lam2', lam=1.)
#         offset2 = pm.Uniform('offset2',upper=float(max(values)))
#         pois1 = pm.Poisson.dist(mu=lam1+offset1)
#         pois2 = pm.Poisson.dist(mu=lam2+offset2)
#         w = pm.Dirichlet('w',a=np.array([1.,1.]))
#         like = pm.Mixture('like',w=w,comp_dists=[pois1,pois2],observed=values)
#     with model:
#         trace = pm.sample(50000,n_init=10000,tune=10000)
#     of1 = np.mean(trace["offset1"][-4000:])
#     of2 = np.mean(trace["offset2"][-4000:])
#     lam1_avg = np.mean(trace['lam1'][-4000:])
#     lam2_avg = np.mean(trace['lam2'][-4000:])
#
#     cluster1Endpoint = int(round(of1 + lam1_avg))
#     cluster2Endpoint = int(round(of2 + lam2_avg))
#
#     return [cluster1Endpoint,cluster2Endpoint]
#
# def model3(values):
#     return [0 for x in range(3)]
#     with pm.Model() as model:
#         lam1 = pm.Exponential('lam1', lam=1.)
#         offset1 = pm.Uniform('offset1',upper=float(max(values)))
#         lam2 = pm.Exponential('lam2', lam=1.)
#         offset2 = pm.Uniform('offset2',upper=float(max(values)))
#         lam3 = pm.Exponential('lam3', lam=1.)
#         offset3 = pm.Uniform('offset3',upper=float(max(values)))
#
#         pois1 = pm.Poisson.dist(mu=lam1+offset1)
#         pois2 = pm.Poisson.dist(mu=lam2+offset2)
#         pois3 = pm.Poisson.dist(mu=lam3+offset3)
#         w = pm.Dirichlet('w',a=np.array([1.,1.,1.]))
#         like = pm.Mixture('like',w=w,comp_dists=[pois1,pois2,pois3],observed=values)
#     with model:
#         trace = pm.sample(70000,n_init=10000,tune=12000)
#     of1 = np.mean(trace["offset1"][-4000:])
#     of2 = np.mean(trace["offset2"][-4000:])
#     of3 = np.mean(trace["offset3"][-4000:])
#     lam1_avg = np.mean(trace['lam1'][-4000:])
#     lam2_avg = np.mean(trace['lam2'][-4000:])
#     lam3_avg = np.mean(trace['lam3'][-4000:])
#
#     cluster1Endpoint = int(round(of1 + lam1_avg))
#     cluster2Endpoint = int(round(of2 + lam2_avg))
#     cluster3Endpoint = int(round(of3 + lam3_avg))
#
#     return [cluster1Endpoint,cluster2Endpoint,cluster3Endpoint]
# def model4(values):
#     return [0 for x in range(4)]
#     with pm.Model() as model:
#         lam1 = pm.Exponential('lam1', lam=1.)
#         offset1 = pm.Uniform('offset1',upper=float(max(values)))
#         lam2 = pm.Exponential('lam2', lam=1.)
#         offset2 = pm.Uniform('offset2',upper=float(max(values)))
#         lam3 = pm.Exponential('lam3', lam=1.)
#         offset3 = pm.Uniform('offset3',upper=float(max(values)))
#         lam4 = pm.Exponential('lam4', lam=1.)
#         offset4 = pm.Uniform('offset4',upper=float(max(values)))
#
#         pois1 = pm.Poisson.dist(mu=lam1+offset1)
#         pois2 = pm.Poisson.dist(mu=lam2+offset2)
#         pois3 = pm.Poisson.dist(mu=lam3+offset3)
#         pois4 = pm.Poisson.dist(mu=lam4+offset4)
#         w = pm.Dirichlet('w',a=np.array([1.,1.,1.,1.]))
#         like = pm.Mixture('like',w=w,comp_dists=[pois1,pois2,pois3,pois4],observed=values)
#     with model:
#         trace = pm.sample(90000,n_init=10000,tune=14000)
#     of1 = np.mean(trace["offset1"][-4000:])
#     of2 = np.mean(trace["offset2"][-4000:])
#     of3 = np.mean(trace["offset3"][-4000:])
#     of4 = np.mean(trace["offset4"][-4000:])
#     lam1_avg = np.mean(trace['lam1'][-4000:])
#     lam2_avg = np.mean(trace['lam2'][-4000:])
#     lam3_avg = np.mean(trace['lam3'][-4000:])
#     lam4_avg = np.mean(trace['lam4'][-4000:])
#
#     cluster1Endpoint = int(round(of1 + lam1_avg))
#     cluster2Endpoint = int(round(of2 + lam2_avg))
#     cluster3Endpoint = int(round(of3 + lam3_avg))
#     cluster4Endpoint = int(round(of4 + lam4_avg))
#
#     return [cluster1Endpoint,cluster2Endpoint,cluster3Endpoint,cluster4Endpoint]
# def model5(values):
#     return [0 for x in range(5)]
#     with pm.Model() as model:
#         lam1 = pm.Exponential('lam1', lam=1.)
#         offset1 = pm.Uniform('offset1',upper=float(max(values)))
#         lam2 = pm.Exponential('lam2', lam=1.)
#         offset2 = pm.Uniform('offset2',upper=float(max(values)))
#         lam3 = pm.Exponential('lam3', lam=1.)
#         offset3 = pm.Uniform('offset3',upper=float(max(values)))
#         lam4 = pm.Exponential('lam4', lam=1.)
#         offset4 = pm.Uniform('offset4',upper=float(max(values)))
#         lam5 = pm.Exponential('lam5', lam=1.)
#         offset5 = pm.Uniform('offset5',upper=float(max(values)))
#
#
#         pois1 = pm.Poisson.dist(mu=lam1+offset1)
#         pois2 = pm.Poisson.dist(mu=lam2+offset2)
#         pois3 = pm.Poisson.dist(mu=lam3+offset3)
#         pois4 = pm.Poisson.dist(mu=lam4+offset4)
#         pois5 = pm.Poisson.dist(mu=lam5+offset5)
#         w = pm.Dirichlet('w',a=np.array([1.,1.,1.,1.,1.]))
#         like = pm.Mixture('like',w=w,comp_dists=[pois1,pois2,pois3,pois4,pois5],observed=values)
#     with model:
#         trace = pm.sample(110000,n_init=10000,tune=16000)
#     of1 = np.mean(trace["offset1"][-4000:])
#     of2 = np.mean(trace["offset2"][-4000:])
#     of3 = np.mean(trace["offset3"][-4000:])
#     of4 = np.mean(trace["offset4"][-4000:])
#     of5 = np.mean(trace["offset5"][-4000:])
#     lam1_avg = np.mean(trace['lam1'][-4000:])
#     lam2_avg = np.mean(trace['lam2'][-4000:])
#     lam3_avg = np.mean(trace['lam3'][-4000:])
#     lam4_avg = np.mean(trace['lam4'][-4000:])
#     lam5_avg = np.mean(trace['lam5'][-4000:])
#
#     cluster1Endpoint = int(round(of1 + lam1_avg))
#     cluster2Endpoint = int(round(of2 + lam2_avg))
#     cluster3Endpoint = int(round(of3 + lam3_avg))
#     cluster4Endpoint = int(round(of4 + lam4_avg))
#     cluster5Endpoint = int(round(of5 + lam5_avg))
#
#     return [cluster1Endpoint,cluster2Endpoint,cluster3Endpoint,cluster4Endpoint,cluster5Endpoint]
# def model6(values):
#     return [0 for x in range(6)]
#     with pm.Model() as model:
#         lam1 = pm.Exponential('lam1', lam=1.)
#         offset1 = pm.Uniform('offset1',upper=float(max(values)))
#         lam2 = pm.Exponential('lam2', lam=1.)
#         offset2 = pm.Uniform('offset2',upper=float(max(values)))
#         lam3 = pm.Exponential('lam3', lam=1.)
#         offset3 = pm.Uniform('offset3',upper=float(max(values)))
#         lam4 = pm.Exponential('lam4', lam=1.)
#         offset4 = pm.Uniform('offset4',upper=float(max(values)))
#         lam5 = pm.Exponential('lam5', lam=1.)
#         offset5 = pm.Uniform('offset5',upper=float(max(values)))
#         lam6 = pm.Exponential('lam6', lam=1.)
#         offset6 = pm.Uniform('offset6',upper=float(max(values)))
#
#
#         pois1 = pm.Poisson.dist(mu=lam1+offset1)
#         pois2 = pm.Poisson.dist(mu=lam2+offset2)
#         pois3 = pm.Poisson.dist(mu=lam3+offset3)
#         pois4 = pm.Poisson.dist(mu=lam4+offset4)
#         pois5 = pm.Poisson.dist(mu=lam5+offset5)
#         pois6 = pm.Poisson.dist(mu=lam6+offset6)
#
#         w = pm.Dirichlet('w',a=np.array([1.,1.,1.,1.,1.,1.]))
#         like = pm.Mixture('like',w=w,comp_dists=[pois1,pois2,pois3,pois4,pois5,pois6],observed=values)
#     with model:
#         trace = pm.sample(130000,n_init=10000,tune=18000)
#     of1 = np.mean(trace["offset1"][-4000:])
#     of2 = np.mean(trace["offset2"][-4000:])
#     of3 = np.mean(trace["offset3"][-4000:])
#     of4 = np.mean(trace["offset4"][-4000:])
#     of5 = np.mean(trace["offset5"][-4000:])
#     of6 = np.mean(trace["offset6"][-4000:])
#     lam1_avg = np.mean(trace['lam1'][-4000:])
#     lam2_avg = np.mean(trace['lam2'][-4000:])
#     lam3_avg = np.mean(trace['lam3'][-4000:])
#     lam4_avg = np.mean(trace['lam4'][-4000:])
#     lam5_avg = np.mean(trace['lam5'][-4000:])
#     lam6_avg = np.mean(trace['lam6'][-4000:])
#
#     cluster1Endpoint = int(round(of1 + lam1_avg))
#     cluster2Endpoint = int(round(of2 + lam2_avg))
#     cluster3Endpoint = int(round(of3 + lam3_avg))
#     cluster4Endpoint = int(round(of4 + lam4_avg))
#     cluster5Endpoint = int(round(of5 + lam5_avg))
#     cluster6Endpoint = int(round(of6 + lam6_avg))
#
#     return [cluster1Endpoint,cluster2Endpoint,cluster3Endpoint,cluster4Endpoint,cluster5Endpoint,cluster6Endpoint]
# def model7(values):
#     return [0 for x in range(7)]
# def model8(values):
#     return [0 for x in range(8)]


plt.rcParams['figure.figsize'] = 14, 6
SEED = 4662 #random.org 
np.random.seed(SEED)
#print('Running on PyMC3 v{}'.format(pm.__version__))
#"%s\t%s\t%s\t%s\t%s\t%d\t%d\t%d\t%d\t%d\t%s\t%s\t%i\n" %(gene, chrom, cds_strand, cds, seq, cds_start, cds_end, num_reads, num_exons, stop_codon, ",".join(stop_positions), ",".join(read_ids),retained_intron))
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
    # print stop_positions
    # print read_ids
    # if gene == "WBGene00018395":
 #        print line
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
            
        # if stop_codon in stop_position_dict:
        #     for x in stop_positions:
        #         stop_position_dict[stop_codon].append(x)
        #     for x in read_ids:
        #         read_id_dict[stop_codon].append(x)
        # else:
        #     stop_position_dict[stop_codon] = stop_positions
        #     read_id_dict[stop_codon] = read_ids
    else:
        for stop_codon2,introns_past_stop2 in stop_position_dict:
            previous_cluster_count = counter 
            stop_pos = stop_position_dict[(stop_codon2,introns_past_stop2)]
            read_ids2 = read_id_dict[(stop_codon2,introns_past_stop2)]
            if len(stop_pos) >= read_threshold:
                #end_points,end_points2 = model(stop_pos,strand = last_strand,read_threshold = read_threshold)
                end_points = model(stop_pos,strand = last_strand,read_threshold = read_threshold)
                
                #print last_gene
                #print end_points
                cluster_sizes.append(len(end_points))
                # if len(end_points) >= 2:
                    # print "here"
                    # sns.distplot(stop_pos,bins=range(min(stop_pos),max(stop_pos)),kde_kws= {'bw':10})
                    # for end_point in end_points:
                    #     plt.arrow(float(end_point),-0.01,0.,0.0075,color="r",head_width=0.01)
                    # for end_point in end_points2:
                    #     plt.arrow(float(end_point),-0.01,0.,0.0075,color="b",head_width=0.01)
                    # plt.title("%s;%d;%s;%d" %(last_gene,len(end_points),last_strand,len(stop_pos)))
                    # plt.show()
                    # time.sleep(1)
                
                for end_point in end_points:
                    if introns_past_stop2 != ():
                        #print last_gene
                        block_starts, block_sizes = convertIntronsToBedBlocks(stop_codon2,end_point,introns_past_stop2,last_strand)
                        # if last_gene == "WBGene00011765":
                        #     print stop_codon2
                        #     print end_point
                        #     print introns_past_stop2
                        #     print strand
                        # for x in block_starts:
                        #     assert x >= 0
                        # for x in block_sizes:
                        #     assert x > 0
                    else:
                        block_starts = [0]
                        if last_strand == '-':
                            block_sizes = [abs(stop_codon2 - int(end_point))]
                        elif last_strand == '+':
                            block_sizes = [abs(stop_codon2 - int(end_point))]
                    if last_strand == '+':
                        #print "here2"
                        #outfile.write("%s\t%d\t%d\t%s\t%d\t%s\n" %(last_chrom,stop_codon2,end_point,last_gene + "-cluster%d" %(counter),1000,last_strand))
                        outfile.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t0,0,255\t%d\t%s\t%s\n" %(last_chrom,stop_codon2,end_point,last_gene + "-cluster%d" %(counter),1000,last_strand,stop_codon2,stop_codon2,len(block_starts),','.join([str(x) for x in block_sizes]),','.join([str(x) for x in block_starts])))
                        # all_clusters.add(last_gene + "-cluster%d" %(counter))
                    elif last_strand == '-':
                        #print "here3"
                        #outfile.write("%s\t%d\t%d\t%s\t%d\t%s\n" %(last_chrom,end_point,stop_codon2,last_gene + "-cluster%d" %(counter),1000,last_strand))
                        outfile.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t0,0,255\t%d\t%s\t%s\n" %(last_chrom,end_point,stop_codon2,last_gene + "-cluster%d" %(counter),1000,last_strand,end_point,end_point,len(block_starts),','.join([str(x) for x in block_sizes]),','.join([str(x) for x in block_starts])))
                        # all_clusters.add(last_gene + "-cluster%d" %(counter))
                    counter += 1
                    # total_counter += 1
                
                # end_points_count = [0 for x in range(len(end_points))]
                
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
                        #assigned_reads_counter += 1
                        assignment_file.write("%s\t%s\n" %(read_id,last_gene + "-cluster%d" %(min_idx + previous_cluster_count)))
                        # all_clusters2.add(last_gene + "-cluster%d" %(min_idx + previous_cluster_count))
                #         end_points_count[min_idx] += 1
                # for i,count in enumerate(end_points_count):
                #     if count == 0:
                #         print last_gene + "-cluster%d" %(i)
                # if len(end_points) != 0 and assigned_reads_counter == 0:
#                     print end_points
#                     print stop_pos
#                     print last_gene
                    
                # sns.distplot(stop_pos,bins=range(min(stop_pos),max(stop_pos)),kde_kws= {'bw':10})
                # plt.title("%s;%d;%d" %(last_gene,len(end_points),len(stop_pos)))
                # plt.show()
                # time.sleep(1)
                # with pm.Model() as model:
                #     lam1 = pm.Exponential('lam1', lam=1.)
                #     offset1 = pm.Uniform('offset1',upper=float(max(new_stop_pos)))
                #     #offset1 = pm.Normal('offset1',mu=0,sd = 100)
                #     lam2 = pm.Exponential('lam2', lam=1.)
                #     #offset2 = pm.Normal('offset2',mu=0,sd = 100)
                #     offset2 = pm.Uniform('offset2',upper= float(max(new_stop_pos)))
                #     pois1 = pm.Poisson.dist(mu=lam1+offset1)
                #     pois2 = pm.Poisson.dist(mu=lam2+offset2)
                #     w = pm.Dirichlet('w',a=np.array([1.,1.]))
                #     like = pm.Mixture('like',w=w,comp_dists=[pois1,pois2],observed=new_stop_pos)
                #
                #
                # with model:
                #     trace = pm.sample(10000,n_init=10000,tune=10000)
                
#                 with pm.Model() as model:
#                     lam1 = pm.Exponential('lam1', lam=1.)
#                     #offset1 = pm.Categorical('offset1',p=[0.1 for x in range(max(new_stop_pos))])
#                     #offset1 = pm.Normal('offset1',mu=0,sd = 100)
#                     lam2 = pm.Exponential('lam2', lam=1.)
#                     #offset2 = pm.Normal('offset2',mu=0,sd = 100)
#                     #offset2 = pm.Categorical('offset2',p=[0.1 for x in range(max(new_stop_pos))])
#                     pois1 = pm.Poisson.dist(mu=lam1+of1)
#                     pois2 = pm.Poisson.dist(mu=lam2+of2)
#                     w = pm.Dirichlet('w',a=np.array([1.,1.]))
#                     like = pm.Mixture('like',w=w,comp_dists=[pois1,pois2],observed=new_stop_pos)
#                 with model:
#                     trace = pm.sample(50000,n_init=10000,tune=100000)
                # of1 = np.mean(trace["offset1"][-2000:])
                # of2 = np.mean(trace["offset2"][-2000:])
                # lam1_avg = np.mean(trace['lam1'][-2000:])
                # lam2_avg = np.mean(trace['lam2'][-2000:])
                #
                # cluster1Endpoint = int(round(min_stop + of1 + lam1_avg))
                # cluster2Endpoint = int(round(min_stop + of2 + lam2_avg))
                #
                # sampled_array = pm.sample_ppc(trace,samples=500,model=model)
                #
                # sampled = []
                # for x in sampled_array['like']:
                #     for y in x:
                #         sampled.append(y)
                #
                # _ = plt.hist(new_stop_pos,density=True, bins =range(max(new_stop_pos)),histtype='step',label='obs')
                # _ = plt.hist(sampled,density=True, bins =range(max(new_stop_pos)),histtype='step',label='model')
                # plt.show()
                # sns.distplot(stop_pos)
#                 plt.title(last_gene)
#                 plt.show()
#                 plt.clf()
                # if last_gene == "WBGene00015682":
                #     sns.distplot(stop_pos,bins=range(min(stop_pos),max(stop_pos)),kde_kws= {'bw':5})
                #     for end_point in end_points:
                #         plt.arrow(float(end_point),-0.01,0.,0.0075,color="b",head_width=0.01)
                #     plt.arrow(float(stop_codon2),-0.01,0.,0.0075,color='r',head_width=0.01)
                #     plt.title("%s;%d;%d" %(last_gene,len(end_points),len(stop_pos)))
                #     plt.show()
                #     time.sleep(1)
                    # flag = True
                    # break
                #putative_offsets = calculatePutativeOffsets(new_stop_pos)
        # if flag:
        #     break
        last_gene = gene
        last_chrom = chrom
        last_strand = strand
        stop_position_dict = {}
        read_id_dict = {}
        counter = 0
        # stop_position_dict[stop_codon] = stop_positions
        # read_id_dict[stop_codon] = read_ids
        stop_position_dict[(stop_codon,introns_past_stop)] = []
        read_id_dict[(stop_codon,introns_past_stop)] = []
        for x in range(len(stop_positions)):
            read_id = read_ids[x]
            stop_pos = stop_positions[x]
            if read_id not in exclusion_reads:
                if (strand == '+' and stop_pos > stop_codon) or (strand == '-' and stop_pos < stop_codon):
                    stop_position_dict[(stop_codon,introns_past_stop)].append(stop_pos)
                    read_id_dict[(stop_codon,introns_past_stop)].append(read_id)
# print max(cluster_sizes)
#print total_counter
#print len(all_clusters)
#print len(all_clusters2)
# plt.hist(cluster_sizes,bins=range(min(cluster_sizes),max(cluster_sizes)+1))
# plt.show()