#!/usr/bin/env python2
import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

if len(sys.argv) > 1:
    mangone_color = sys.argv[1]
    jan_color = sys.argv[2]
    our_color = sys.argv[3]
else:
    mangone_color = '#B3DE8E'
    jan_color = '#F99B9B'
    our_color = '#A7CEE2'

font = {"family":"sans-serif",
        "sans-serif":["Helvetica"],
        "size":10}
matplotlib.rc('font',**font)

l1_in  = open("../../../results/utrs/beds/L1_utrs.bed")
l2_in  = open("../../../results/utrs/beds/L2_utrs.bed")
l3_in  = open("../../../results/utrs/beds/L3_utrs.bed")
l4_in  = open("../../../results/utrs/beds/L4_utrs.bed")
ya_in  = open("../../../results/utrs/beds/YA_utrs.bed")
ga_in  = open("../../../results/utrs/beds/GA_utrs.bed")
ml_in  = open("../../../results/utrs/beds/ML_utrs.bed")
all_in = open("../../../results/utrs/beds/all_isoforms_utrs.bed")
mangone_in = open("../../../references/utrs/mangone_uniq_utrs.bed")
jan_in     = open("../../../references/utrs/jan_uniq_utrs.bed")


l1_lengths = []
l2_lengths = []
l3_lengths = []
l4_lengths = []
ya_lengths = []
ga_lengths = []
ml_lengths = []
all_lengths = []
mangone_lengths = []
jan_lengths = []
for line in l1_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    l1_lengths.append(end - start)
    
for line in l2_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    l2_lengths.append(end - start)
    
for line in l3_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    l3_lengths.append(end - start)
    
for line in l4_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    l4_lengths.append(end - start)
    
for line in ya_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    ya_lengths.append(end - start)
    
for line in ga_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    ga_lengths.append(end - start)
    
for line in ml_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    ml_lengths.append(end - start)
    
for line in all_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    all_lengths.append(end - start)
    
for line in mangone_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    mangone_lengths.append(end - start)
for line in jan_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    jan_lengths.append(end - start)
    
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 16}
    
matplotlib.rc("font",**font)

sns.set(style="white",font="Helvetica")
plt.figure(num=None,figsize=(4.166,2.5))
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
sns.kdeplot(jan_lengths,label = "Jan et al 3'UTRs",gridsize=10000,color=jan_color)
sns.kdeplot(mangone_lengths,label = "Mangone et al 3'UTRs",gridsize=10000,color=mangone_color)
sns.kdeplot(all_lengths,label = "This study",gridsize=10000,color=our_color)

plt.legend(fontsize=10,frameon=False)

plt.ylim(0,0.006)
plt.xlim(0,1000)
plt.xticks(fontsize=10)
plt.yticks(fontsize=10)
plt.xlabel("3'UTR Length",fontsize=10)
plt.ylabel("Proportion of 3'UTRs",fontsize=10)
plt.tight_layout()
plt.savefig("../../../figures/figure3/figure3D.pdf")
plt.clf()

