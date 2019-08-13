#!/usr/bin/env python

import sys
import pysam
import matplotlib
import matplotlib.pyplot as plt

font = {"family":"sans-serif",
        "sans-serif":["Helvetica"],
        "size":8}
matplotlib.rc('font',**font)

infile1 = pysam.AlignmentFile("../../data/L1/combined/L1.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/L1/combined/L1_beta_filtered.bam",'rb')

lengths1 = []
for read in infile1.fetch():
    lengths1.append(float(len(read.query_sequence))/1000.)

lengths2 = []
for read in infile2.fetch():
    lengths2.append(float(len(read.query_sequence))/1000.)


plt.figure(num=None,figsize=(3,1.6))

# plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All L1 Reads")
# plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length L1 Reads")
plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All L1 Reads")
plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length L1 Reads")

plt.legend()
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# plt.show()
plt.ylabel("Number of Reads")
plt.tight_layout()
plt.savefig("../../figures/supplementals/sfigure3/L1_read_length_distribution.png", dpi=300)


infile1 = pysam.AlignmentFile("../../data/L2/combined/L2.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/L2/combined/L2_beta_filtered.bam",'rb')

lengths1 = []
for read in infile1.fetch():
    lengths1.append(float(len(read.query_sequence))/1000.)

lengths2 = []
for read in infile2.fetch():
    lengths2.append(float(len(read.query_sequence))/1000.)


plt.figure(num=None,figsize=(3,1.6))

plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All L2 Reads")
plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length L2 Reads")
plt.legend()
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# plt.show()
plt.ylabel("Number of Reads")
plt.tight_layout()
plt.savefig("../../figures/supplementals/sfigure3/L2_read_length_distribution.png", dpi=300)



infile1 = pysam.AlignmentFile("../../data/L3/combined/L3.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/L3/combined/L3_beta_filtered.bam",'rb')

lengths1 = []
for read in infile1.fetch():
    lengths1.append(float(len(read.query_sequence))/1000.)

lengths2 = []
for read in infile2.fetch():
    lengths2.append(float(len(read.query_sequence))/1000.)


plt.figure(num=None,figsize=(3,1.6))

plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All L3 Reads")
plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length L3 Reads")
plt.legend()
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# plt.show()
plt.ylabel("Number of Reads")
plt.tight_layout()
plt.savefig("../../figures/supplementals/sfigure3/L3_read_length_distribution.png", dpi=300)


infile1 = pysam.AlignmentFile("../../data/L4/combined/L4.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/L4/combined/L4_beta_filtered.bam",'rb')

lengths1 = []
for read in infile1.fetch():
    lengths1.append(float(len(read.query_sequence))/1000.)

lengths2 = []
for read in infile2.fetch():
    lengths2.append(float(len(read.query_sequence))/1000.)


plt.figure(num=None,figsize=(3,1.6))

plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All L4 Reads")
plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length L4 Reads")
plt.legend()
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# plt.show()
plt.ylabel("Number of Reads")
plt.tight_layout()
plt.savefig("../../figures/supplementals/sfigure3/L4_read_length_distribution.png", dpi=300)


infile1 = pysam.AlignmentFile("../../data/young_adult/combined/YA.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/young_adult/combined/YA_beta_filtered.bam",'rb')

lengths1 = []
for read in infile1.fetch():
    lengths1.append(float(len(read.query_sequence))/1000.)

lengths2 = []
for read in infile2.fetch():
    lengths2.append(float(len(read.query_sequence))/1000.)


plt.figure(num=None,figsize=(3,1.6))

plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All yAd Reads")
plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length yAd Reads")
plt.legend()
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# plt.show()
plt.ylabel("Number of Reads")
plt.tight_layout()
plt.savefig("../../figures/supplementals/sfigure3/YA_read_length_distribution.png", dpi=300)

infile1 = pysam.AlignmentFile("../../data/adult/combined/GA.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/adult/combined/GA_beta_filtered.bam",'rb')

lengths1 = []
for read in infile1.fetch():
    lengths1.append(float(len(read.query_sequence))/1000.)

lengths2 = []
for read in infile2.fetch():
    lengths2.append(float(len(read.query_sequence))/1000.)


plt.figure(num=None,figsize=(3,1.6))

plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All mAd Reads")
plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length mAd Reads")
plt.legend()
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# plt.show()
plt.ylabel("Number of Reads")
plt.tight_layout()
# plt.xticks([])
plt.savefig("../../figures/supplementals/sfigure3/GA_read_length_distribution.png", dpi=300)

infile1 = pysam.AlignmentFile("../../data/male/combined/ML.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/male/combined/ML_beta_filtered.bam",'rb')

lengths1 = []
for read in infile1.fetch():
    lengths1.append(float(len(read.query_sequence))/1000.)

lengths2 = []
for read in infile2.fetch():
    lengths2.append(float(len(read.query_sequence))/1000.)


plt.figure(num=None,figsize=(3,1.6))

plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All male Reads")
plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length male Reads")
plt.legend()
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# plt.show()
plt.ylabel("Number of Reads")
# plt.xlabel("Read length (nt)")
plt.tight_layout()
plt.savefig("../../figures/supplementals/sfigure3/ML_read_length_distribution.png", dpi=300)