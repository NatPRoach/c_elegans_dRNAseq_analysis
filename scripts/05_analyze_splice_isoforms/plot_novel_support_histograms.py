#!/usr/bin/env python2

import sys
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

def getSupportCounts(infile):
    read_counts =[]
    for line in infile:
        fields = line.strip().split()
        read_counts.append(np.log2(float(fields[7])))
    count = 0 
    for x in read_counts:
        if x >= 2.:
            count += 1
    print float(count) * 100 / len(read_counts)
    return read_counts


if len(sys.argv) > 1:
    color = sys.argv[1]
else:
    color = "#6A4098"
font = {"family":"sans-serif",
        "sans-serif":["Helvetica"],
        "size":10}
matplotlib.rc('font',**font)

### Uncomment out below to get support histograms broken up by stage.
# l1_counts = getSupportCounts(open("../../results/new_isoforms/L1_new_isoforms.tsv"))
# plt.figure(num=None,figsize=(4.1667,2.5))
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.hist(l1_counts,bins= [x - 0.5 for x in range(int(np.round(min(l1_counts))),8)],density=True,rwidth=0.4,color=color,edgecolor="black")
# plt.title("L1")
# plt.xlabel("Log2 Read Counts")
# plt.ylabel("proportion of novel splice isoforms")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/novelSupportHistograms/01_l1NovelIsoformSupport.pdf")
# plt.clf()
#
# l2_counts = getSupportCounts(open("../../results/new_isoforms/L2_new_isoforms.tsv"))
# plt.figure(num=None,figsize=(4.1667,2.5))
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.hist(l2_counts,bins= [x - 0.5 for x in range(int(np.round(min(l2_counts))),8)],density=True,rwidth=0.4,color=color,edgecolor="black")
# plt.title("L2")
# plt.xlabel("Log2 Read Counts")
# plt.ylabel("proportion of novel splice isoforms")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/novelSupportHistograms/02_l2NovelIsoformSupport.pdf")
# plt.clf()
#
# l3_counts = getSupportCounts(open("../../results/new_isoforms/L3_new_isoforms.tsv"))
# plt.figure(num=None,figsize=(4.1667,2.5))
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.hist(l3_counts,bins= [x - 0.5 for x in range(int(np.round(min(l3_counts))),8)],density=True,rwidth=0.4,color=color,edgecolor="black")
# plt.title("L3")
# plt.xlabel("Log2 Read Counts")
# plt.ylabel("proportion of novel splice isoforms")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/novelSupportHistograms/03_l3NovelIsoformSupport.pdf")
# plt.clf()
#
# l4_counts = getSupportCounts(open("../../results/new_isoforms/L4_new_isoforms.tsv"))
# plt.figure(num=None,figsize=(4.1667,2.5))
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.hist(l4_counts,bins= [x - 0.5 for x in range(int(np.round(min(l4_counts))),8)],density=True,rwidth=0.4,color=color,edgecolor="black")
# plt.title("L4")
# plt.xlabel("Log2 Read Counts")
# plt.ylabel("proportion of novel splice isoforms")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/novelSupportHistograms/04_l4NovelIsoformSupport.pdf")
# plt.clf()
#
# ya_counts = getSupportCounts(open("../../results/new_isoforms/YA_new_isoforms.tsv"))
# plt.figure(num=None,figsize=(4.1667,2.5))
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.hist(ya_counts,bins= [x - 0.5 for x in range(int(np.round(min(ya_counts))),8)],density=True,rwidth=0.4,color=color,edgecolor="black")
# plt.title("young adult")
# plt.xlabel("Log2 Read Counts")
# plt.ylabel("proportion of novel splice isoforms")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/novelSupportHistograms/05_yaNovelIsoformSupport.pdf")
# plt.clf()
#
# ga_counts = getSupportCounts(open("../../results/new_isoforms/GA_new_isoforms.tsv"))
# plt.figure(num=None,figsize=(4.1667,2.5))
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.hist(ga_counts,bins= [x - 0.5 for x in range(int(np.round(min(ga_counts))),8)],density=True,rwidth=0.4,color=color,edgecolor="black")
# plt.title("mature adult")
# plt.xlabel("Log2 Read Counts")
# plt.ylabel("proportion of novel splice isoforms")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/novelSupportHistograms/06_gaNovelIsoformSupport.pdf")
# plt.clf()
#
# ml_counts = getSupportCounts(open("../../results/new_isoforms/ML_new_isoforms.tsv"))
# plt.figure(num=None,figsize=(4.1667,2.5))
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.hist(ml_counts,bins= [x - 0.5 for x in range(int(np.round(min(ml_counts))),8)],density=True,rwidth=0.4,color=color,edgecolor="black")
# plt.title("male")
# plt.xlabel("Log2 Read Counts")
# plt.ylabel("proportion of novel splice isoforms")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/novelSupportHistograms/07_mlNovelIsoformSupport.pdf")
# plt.clf()
#
# all_counts = getSupportCounts(open("../../results/new_isoforms/all_new_isoforms.tsv"))
# plt.figure(num=None,figsize=(4.1667,2.5))
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.hist(all_counts,bins= [x - 0.5 for x in range(int(np.round(min(all_counts))),8)],density=True,rwidth=0.4,color=color,edgecolor="black")
# plt.title("all")
# plt.xlabel("Log2 Read Counts")
# plt.ylabel("proportion of novel splice isoforms")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/novelSupportHistograms/08_allNovelIsoformSupport.pdf")
# plt.clf()

# all_counts = getSupportCounts(open("../../results/new_isoforms/all_new_isoforms.tsv"))
# plt.figure(num=None,figsize=(4.1667,2.5))
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.hist(all_counts,bins= [x - 0.5 for x in range(int(np.round(min(all_counts))),8)],density=True,rwidth=0.4,color=color,edgecolor="black")
# plt.xlabel("Log2 Read Counts")
# plt.ylabel("proportion of\nnovel splice isoforms")
# plt.ylim((0,0.6))
# plt.tight_layout()
# plt.savefig("../../figures/figure2/figure2F.pdf")
# plt.clf()

# ### Sensitive
# all_counts = getSupportCounts(open("../../results/new_isoforms/all_sensitive_new_isoforms.tsv"))
# plt.figure(num=None,figsize=(4.1667,2.5))
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.hist(all_counts,bins= [x - 0.5 for x in range(int(np.round(min(all_counts))),8)],density=True,rwidth=0.4,color=color,edgecolor="black")
# plt.xticks([0,1,2,3,4,5,6])
# plt.xlabel("Log2 Read Counts")
# plt.ylabel("proportion of\nnovel splice isoforms")
# plt.ylim((0,0.6))
# plt.tight_layout()
# plt.savefig("../../figures/figure2/sensitive_figure2F.pdf")
# plt.clf()

### Stringent
all_counts = getSupportCounts(open("../../results/new_isoforms/all_stringent_new_isoforms.tsv"))
plt.figure(num=None,figsize=(4.1667,2.5))
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.hist(all_counts,bins= [x - 0.5 for x in range(int(np.round(min(all_counts))),8)],density=True,rwidth=0.4,color=color,edgecolor="black")
plt.xticks([0,1,2,3,4,5,6])
plt.xlabel("Log2 Read Counts")
plt.ylabel("proportion of\nnovel splice isoforms")
plt.ylim((0,0.6))
plt.tight_layout()
plt.savefig("../../figures/figure2/figure2F.pdf")
plt.clf()