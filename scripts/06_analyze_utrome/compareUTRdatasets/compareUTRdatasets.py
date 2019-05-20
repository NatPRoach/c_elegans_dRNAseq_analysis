#!/usr/bin/env python2
## Script written by Nathan Roach ( nroach2@jhu.edu )
## Goal is to compare our dataset with two previously published UTR datasets
## Mangone et al, and Jan ... Bartel et al 
## Actually something of a tricky problem, need to compare the sets, with some degree of tolerance, to allow for slight differences
## Currently dealing with this by comparing everything vs a master set



import sets

my_utrs = open("../../../results/utrs/beds/all_isoforms_utrs.bed",'r')
mangone_utrs  = open("../../../references/utrs/mangone_utrs.bed",'r')
bartel_utrs   = open("../../../references/utrs/jan_utrs.bed",'r')
wormbase_utrs = open("../../../references/utrs/wormbase_utrs.bed",'r')
outfile = open("../../../results/overlaps/utr.overlap.matrix",'w')
outfile2 = open("../../../results/utrs/beds/novel/all_novel_utrs.bed",'w')
window = 10

my_utr_dict = {}
recall_dict = {}
my_bed_dict = {}
for line in my_utrs:
    fields = line.strip().split()
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    cluster_id = fields[3]
    strand = fields[5]
    if strand == '+':
        id_tuple = (chrom,strand,start)
        id_tuple2 = (chrom,strand,start,end)
        cleavage_site = end
    elif strand == '-':
        id_tuple = (chrom,strand,end)
        id_tuple2 = (chrom,strand,end,start)
        cleavage_site = start
    else:
        print "Something went wrong"
    if id_tuple in my_utr_dict:
        my_utr_dict[id_tuple].append(cleavage_site)
    else:
        my_utr_dict[id_tuple] = [cleavage_site]
    if id_tuple2 in recall_dict:
        recall_dict[id_tuple2].append(cluster_id)
    else:
        recall_dict[id_tuple2] = [cluster_id]
    my_bed_dict[cluster_id] = line
    
mangone_utr_dict = {}
for line in mangone_utrs:
    fields = line.strip().split()
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    strand = fields[5]
    if strand == '+':
        id_tuple = (chrom,strand,start)
        cleavage_site = end
    elif strand == '-':
        id_tuple = (chrom,strand,end)
        cleavage_site = start
    else:
        print "Something went wrong"
    if id_tuple in mangone_utr_dict:
        mangone_utr_dict[id_tuple].append(cleavage_site)
    else:
        mangone_utr_dict[id_tuple] = [cleavage_site]
    
bartel_utr_dict = {}
for line in bartel_utrs:
    fields = line.strip().split()
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    strand = fields[5]
    if strand == '+':
        id_tuple = (chrom,strand,start)
        cleavage_site = end
    elif strand == '-':
        id_tuple = (chrom,strand,end)
        cleavage_site = start
    else:
        print "Something went wrong"
    if id_tuple in bartel_utr_dict:
        bartel_utr_dict[id_tuple].append(cleavage_site)
    else:
        bartel_utr_dict[id_tuple] = [cleavage_site]


wormbase_utr_dict = {}
for line in wormbase_utrs:
    fields = line.strip().split()
    chrom = fields[0]
    start = int(fields[1])
    end = int(fields[2])
    strand = fields[5]
    if strand == '+':
        id_tuple = (chrom,strand,start)
        cleavage_site = end
    elif strand == '-':
        id_tuple = (chrom,strand,end)
        cleavage_site = start
    else:
        print "Something went wrong"
    if id_tuple in wormbase_utr_dict:
        wormbase_utr_dict[id_tuple].append(cleavage_site)
    else:
        wormbase_utr_dict[id_tuple] = [cleavage_site]

total_set = sets.Set()
my_set = sets.Set()
mangone_set = sets.Set()
bartel_set = sets.Set()
# wormbase_set = sets.Set()

for chrom,strand,stop_codon in my_utr_dict:
    for cleavage_site in my_utr_dict[(chrom,strand,stop_codon)]:
        #if (chrom,strand,stop_codon,cleavage_site) in total_set:
        #    print chrom,strand,stop_codon,cleavage_site
        total_set.add((chrom,strand,stop_codon,cleavage_site))
        my_set.add((chrom,strand,stop_codon,cleavage_site))

for chrom,strand,stop_codon in mangone_utr_dict:
    for cleavage_site in mangone_utr_dict[(chrom,strand,stop_codon)]:
        if (chrom,strand,stop_codon,cleavage_site) in total_set:
            mangone_set.add((chrom,strand,stop_codon,cleavage_site))
        else:
            flag = True
            for x in range(cleavage_site - window, cleavage_site + window + 1): # is it within a window? 
                if (chrom,strand,stop_codon,x) in total_set:
                    mangone_set.add((chrom,strand,stop_codon,x))
                    flag = False
                    break
            if flag:
                total_set.add((chrom,strand,stop_codon,cleavage_site))
                mangone_set.add((chrom,strand,stop_codon,cleavage_site))


for chrom,strand,stop_codon in bartel_utr_dict:
    for cleavage_site in bartel_utr_dict[(chrom,strand,stop_codon)]:
        if (chrom,strand,stop_codon,cleavage_site) in total_set:
            bartel_set.add((chrom,strand,stop_codon,cleavage_site))
        else:
            flag = True
            for x in range(cleavage_site - window, cleavage_site + window + 1): # is it within a window? 
                if (chrom,strand,stop_codon,x) in total_set:
                    bartel_set.add((chrom,strand,stop_codon,x))
                    flag = False
                    break
            if flag:
                total_set.add((chrom,strand,stop_codon,cleavage_site))
                bartel_set.add((chrom,strand,stop_codon,cleavage_site))


# for chrom,strand,stop_codon in wormbase_utr_dict:
#     for cleavage_site in wormbase_utr_dict[(chrom,strand,stop_codon)]:
#         if (chrom,strand,stop_codon,cleavage_site) in total_set:
#             wormbase_set.add((chrom,strand,stop_codon,cleavage_site))
#         else:
#             flag = True
#             for x in range(cleavage_site - window, cleavage_site + window + 1): # is it within a window?
#                 # if stop_codon == 11420742:
#                 #     print "here1" , x
#                 #     if x == 11420485:
#                 #         print "here"
#                 if (chrom,strand,stop_codon,x) in total_set:
#                     wormbase_set.add((chrom,strand,stop_codon,x))
#                     flag = False
#                     break
#             if flag:
#                 total_set.add((chrom,strand,stop_codon,cleavage_site))
#                 wormbase_set.add((chrom,strand,stop_codon,cleavage_site))

for t in total_set:
    # outfile.write("%i\t%i\t%i\t%i\n" %(t in my_set, t in mangone_set,t in bartel_set,t in wormbase_set))
    outfile.write("%i\t%i\t%i\n" %(t in my_set, t in mangone_set,t in bartel_set))
    # if t in my_set and t not in mangone_set  and t not in bartel_set and t not in wormbase_set:
    #     for cluster_id in recall_dict[t]:
    #         outfile2.write(my_bed_dict[cluster_id])

total_set2 = sets.Set()
for chrom,strand,stop_codon in mangone_utr_dict:
    for cleavage_site in mangone_utr_dict[(chrom,strand,stop_codon)]:
        total_set2.add((chrom,strand,stop_codon,cleavage_site))
for chrom,strand,stop_codon in bartel_utr_dict:
    for cleavage_site in bartel_utr_dict[(chrom,strand,stop_codon)]:
        total_set2.add((chrom,strand,stop_codon,cleavage_site))
for chrom,strand,stop_codon in wormbase_utr_dict:
    for cleavage_site in wormbase_utr_dict[(chrom,strand,stop_codon)]:
        total_set2.add((chrom,strand,stop_codon,cleavage_site))

for chrom,strand,stop_codon in my_utr_dict:
    for cleavage_site in my_utr_dict[(chrom,strand,stop_codon)]:
        if (chrom,strand,stop_codon,cleavage_site) in total_set2:
            continue
            #wormbase_set.add((chrom,strand,stop_codon,cleavage_site))
        else:
            flag = True
            for x in range(cleavage_site - window, cleavage_site + window + 1): # is it within a window?
                # if stop_codon == 11420742:
                #     print "here1" , x
                #     if x == 11420485:
                #         print "here"
                if (chrom,strand,stop_codon,x) in total_set2:
                    flag = False
                    break
            if flag:
                t = (chrom,strand,stop_codon,cleavage_site)
                for cluster_id in recall_dict[t]:
                    outfile2.write(my_bed_dict[cluster_id])
