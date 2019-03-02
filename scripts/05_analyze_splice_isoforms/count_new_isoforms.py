#!/usr/bin/env python2


import sys
import sets

l1_novel_in = open("results/realigned_new_isoforms/L1_new_isoforms.tsv")
l2_novel_in = open("results/realigned_new_isoforms/L2_new_isoforms.tsv")
l3_novel_in = open("results/realigned_new_isoforms/L3_new_isoforms.tsv")
l4_novel_in = open("results/realigned_new_isoforms/L4_new_isoforms.tsv")
ya_novel_in = open("results/realigned_new_isoforms/YA_new_isoforms.tsv")
ga_novel_in = open("results/realigned_new_isoforms/GA_new_isoforms.tsv")
ml_novel_in = open("results/realigned_new_isoforms/ML_new_isoforms.tsv")

l1_novel_set = sets.Set()
l2_novel_set = sets.Set()
l3_novel_set = sets.Set()
l4_novel_set = sets.Set()
ya_novel_set = sets.Set()
ga_novel_set = sets.Set()
ml_novel_set = sets.Set()

for line in l1_novel_in:
    fields = line.strip().split()
    if len(fields) == 16:
        if fields[14] != '1':
            intron_string = fields[15]
            introns = [fields[1]]
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            l1_novel_set.add(introns)

for line in l2_novel_in:
    fields = line.strip().split()
    if len(fields) == 16:
        if fields[14] != '1':
            intron_string = fields[15]
            introns = [fields[1]]
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            l2_novel_set.add(introns)
            
for line in l3_novel_in:
    fields = line.strip().split()
    if len(fields) == 16:
        if fields[14] != '1':
            intron_string = fields[15]
            introns = [fields[1]]
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            l3_novel_set.add(introns)
            
for line in l4_novel_in:
    fields = line.strip().split()
    if len(fields) == 16:
        if fields[14] != '1':
            intron_string = fields[15]
            introns = [fields[1]]
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            l4_novel_set.add(introns)
            
for line in ya_novel_in:
    fields = line.strip().split()
    if len(fields) == 16:
        if fields[14] != '1':
            intron_string = fields[15]
            introns = [fields[1]]
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            ya_novel_set.add(introns)

for line in ga_novel_in:
    fields = line.strip().split()
    if len(fields) == 16:
        if fields[14] != '1':
            intron_string = fields[15]
            introns = [fields[1]]
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            ga_novel_set.add(introns)
            
for line in ml_novel_in:
    fields = line.strip().split()
    if len(fields) == 16:
        if fields[14] != '1':
            intron_string = fields[15]
            introns = [fields[1]]
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            ml_novel_set.add(introns)


all_novel_set = l1_novel_set.union(l2_novel_set.union(l3_novel_set.union(l4_novel_set.union(ya_novel_set.union(ga_novel_set.union(ml_novel_set))))))
print "Novel Intron Chain Isoforms"
print "L1, %d" %(len(l1_novel_set))
print "L2, %d" %(len(l2_novel_set))
print "L3, %d" %(len(l3_novel_set))
print "L4, %d" %(len(l4_novel_set))
print "YA, %d" %(len(ya_novel_set))
print "GA, %d" %(len(ga_novel_set))
print "ML, %d" %(len(ml_novel_set))
print "All, %d" %(len(all_novel_set))
