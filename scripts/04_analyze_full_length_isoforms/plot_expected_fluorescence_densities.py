#!/usr/bin/env python

import sys
import pysam
import matplotlib
import matplotlib.pyplot as plt
import numpy as np

font = {"family":"sans-serif",
        "sans-serif":["Helvetica"],
        "size":8}
matplotlib.rc('font',**font)

def plotExpectedFluorescenceDensities(infile1,infile2,label,outfile):
    lengths1 = []
    weights1 = []
    for read in infile1.fetch():
        # for _ in range(len(read.query_sequence)):
        lengths1.append(float(len(read.query_sequence))/1000.)
        weights1.append(len(read.query_sequence))

    lengths2 = []
    weights2 = []
    for read in infile2.fetch():
        # for _ in range(len(read.query_sequence)):
        lengths2.append(float(len(read.query_sequence))/1000.)
        weights2.append(len(read.query_sequence))

    plt.figure(num=None,figsize=(3,1.6))

    # plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All L1 Reads")
    # plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length L1 Reads")
    plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],weights=weights1,alpha=0.5,label="All " + label,density=True)
    plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],weights=weights2,alpha=0.4,label="Full-Length\n" + label,density=True)

    plt.legend(loc="upper right")
    ax = plt.gca()
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    plt.xscale("log",basex=2)
    # plt.xticks([np.log2(0.025),np.log2(0.2),np.log2(0.5),np.log2(1),np.log2(2),np.log2(4),np.log2(6)],[".025",".2",".5","1","2","4","6"])
    plt.xticks([0.025,0.2,0.5,1,2,4,6],[".025",".2",".5","1","2","4","6"])
    plt.xlim([0.01,120])
    plt.ylim([0,1.5])
    plt.yticks([0,0.5,1.0])
# plt.show()
    plt.ylabel("Expected fluorescence density")
    plt.tight_layout()
    plt.savefig(outfile, dpi=300)

infile1 = pysam.AlignmentFile("../../data/L1/combined/L1.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/L1/combined/L1_beta_filtered.bam",'rb')
plotExpectedFluorescenceDensities(infile1,infile2,"L1 Reads","../../figures/supplementals/sfigure3/L1_expected_fluorescence_distribution.png")

infile1 = pysam.AlignmentFile("../../data/L1/combined/L1.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/L1/combined/L1_beta_filtered.bam",'rb')
plotExpectedFluorescenceDensities(infile1,infile2,"L1 Reads","../../figures/supplementals/sfigure3/L1_expected_fluorescence_distribution.png")


# lengths1 = []
# for read in infile1.fetch():
#     lengths1.append(float(len(read.query_sequence))/1000.)
#
# lengths2 = []
# for read in infile2.fetch():
#     lengths2.append(float(len(read.query_sequence))/1000.)
#
#
# plt.figure(num=None,figsize=(3,1.6))
#
# # plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All L1 Reads")
# # plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length L1 Reads")
# plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All L1 Reads")
# plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length L1 Reads")
#
# plt.legend()
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# # plt.show()
# plt.ylabel("Number of Reads")
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/sfigure3/L1_read_length_distribution.png", dpi=300)

infile1 = pysam.AlignmentFile("../../data/L2/combined/L2.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/L2/combined/L2_beta_filtered.bam",'rb')
plotExpectedFluorescenceDensities(infile1,infile2,"L2 Reads","../../figures/supplementals/sfigure3/L2_expected_fluorescence_distribution.png")

# infile1 = pysam.AlignmentFile("../../data/L2/combined/L2.original_alignment.bam",'rb')
# infile2 = pysam.AlignmentFile("../../data/L2/combined/L2_beta_filtered.bam",'rb')
#
# lengths1 = []
# for read in infile1.fetch():
#     lengths1.append(float(len(read.query_sequence))/1000.)
#
# lengths2 = []
# for read in infile2.fetch():
#     lengths2.append(float(len(read.query_sequence))/1000.)
#
#
# plt.figure(num=None,figsize=(3,1.6))
#
# plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All L2 Reads")
# plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length L2 Reads")
# plt.legend()
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# # plt.show()
# plt.ylabel("Number of Reads")
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/sfigure3/L2_read_length_distribution.png", dpi=300)
#
#

infile1 = pysam.AlignmentFile("../../data/L3/combined/L3.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/L3/combined/L3_beta_filtered.bam",'rb')
plotExpectedFluorescenceDensities(infile1,infile2,"L3 Reads","../../figures/supplementals/sfigure3/L3_expected_fluorescence_distribution.png")

# infile1 = pysam.AlignmentFile("../../data/L3/combined/L3.original_alignment.bam",'rb')
# infile2 = pysam.AlignmentFile("../../data/L3/combined/L3_beta_filtered.bam",'rb')
#
# lengths1 = []
# for read in infile1.fetch():
#     lengths1.append(float(len(read.query_sequence))/1000.)
#
# lengths2 = []
# for read in infile2.fetch():
#     lengths2.append(float(len(read.query_sequence))/1000.)
#
#
# plt.figure(num=None,figsize=(3,1.6))
#
# plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All L3 Reads")
# plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length L3 Reads")
# plt.legend()
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# # plt.show()
# plt.ylabel("Number of Reads")
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/sfigure3/L3_read_length_distribution.png", dpi=300)
#

infile1 = pysam.AlignmentFile("../../data/L4/combined/L4.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/L4/combined/L4_beta_filtered.bam",'rb')
plotExpectedFluorescenceDensities(infile1,infile2,"L4 Reads","../../figures/supplementals/sfigure3/L4_expected_fluorescence_distribution.png")

# infile1 = pysam.AlignmentFile("../../data/L4/combined/L4.original_alignment.bam",'rb')
# infile2 = pysam.AlignmentFile("../../data/L4/combined/L4_beta_filtered.bam",'rb')
#
# lengths1 = []
# for read in infile1.fetch():
#     lengths1.append(float(len(read.query_sequence))/1000.)
#
# lengths2 = []
# for read in infile2.fetch():
#     lengths2.append(float(len(read.query_sequence))/1000.)
#
#
# plt.figure(num=None,figsize=(3,1.6))
#
# plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All L4 Reads")
# plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length L4 Reads")
# plt.legend()
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# # plt.show()
# plt.ylabel("Number of Reads")
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/sfigure3/L4_read_length_distribution.png", dpi=300)
#

infile1 = pysam.AlignmentFile("../../data/young_adult/combined/YA.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/young_adult/combined/YA_beta_filtered.bam",'rb')

plotExpectedFluorescenceDensities(infile1,infile2,"yAd Reads","../../figures/supplementals/sfigure3/YA_expected_fluorescence_distribution.png")

# infile1 = pysam.AlignmentFile("../../data/young_adult/combined/YA.original_alignment.bam",'rb')
# infile2 = pysam.AlignmentFile("../../data/young_adult/combined/YA_beta_filtered.bam",'rb')
#
# lengths1 = []
# for read in infile1.fetch():
#     lengths1.append(float(len(read.query_sequence))/1000.)
#
# lengths2 = []
# for read in infile2.fetch():
#     lengths2.append(float(len(read.query_sequence))/1000.)
#
#
# plt.figure(num=None,figsize=(3,1.6))
#
# plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All yAd Reads")
# plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length yAd Reads")
# plt.legend()
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# # plt.show()
# plt.ylabel("Number of Reads")
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/sfigure3/YA_read_length_distribution.png", dpi=300)

infile1 = pysam.AlignmentFile("../../data/adult/combined/GA.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/adult/combined/GA_beta_filtered.bam",'rb')
plotExpectedFluorescenceDensities(infile1,infile2,"mAd Reads","../../figures/supplementals/sfigure3/GA_expected_fluorescence_distribution.png")

# infile1 = pysam.AlignmentFile("../../data/adult/combined/GA.original_alignment.bam",'rb')
# infile2 = pysam.AlignmentFile("../../data/adult/combined/GA_beta_filtered.bam",'rb')
#
# lengths1 = []
# for read in infile1.fetch():
#     lengths1.append(float(len(read.query_sequence))/1000.)
#
# lengths2 = []
# for read in infile2.fetch():
#     lengths2.append(float(len(read.query_sequence))/1000.)
#
#
# plt.figure(num=None,figsize=(3,1.6))
#
# plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All mAd Reads")
# plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length mAd Reads")
# plt.legend()
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# # plt.show()
# plt.ylabel("Number of Reads")
# plt.tight_layout()
# # plt.xticks([])
# plt.savefig("../../figures/supplementals/sfigure3/GA_read_length_distribution.png", dpi=300)

infile1 = pysam.AlignmentFile("../../data/male/combined/ML.original_alignment.bam",'rb')
infile2 = pysam.AlignmentFile("../../data/male/combined/ML_beta_filtered.bam",'rb')
plotExpectedFluorescenceDensities(infile1,infile2,"male Reads","../../figures/supplementals/sfigure3/ML_expected_fluorescence_distribution.png")

# infile1 = pysam.AlignmentFile("../../data/male/combined/ML.original_alignment.bam",'rb')
# infile2 = pysam.AlignmentFile("../../data/male/combined/ML_beta_filtered.bam",'rb')
#
# lengths1 = []
# for read in infile1.fetch():
#     lengths1.append(float(len(read.query_sequence))/1000.)
#
# lengths2 = []
# for read in infile2.fetch():
#     lengths2.append(float(len(read.query_sequence))/1000.)
#
#
# plt.figure(num=None,figsize=(3,1.6))
#
# plt.hist(lengths1,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.5,label="All male Reads")
# plt.hist(lengths2,bins=[float(x - 0.5)/1000. for x in range(10, 6000)],alpha=0.4,label="Full-Length male Reads")
# plt.legend()
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# # plt.show()
# plt.ylabel("Number of Reads")
# # plt.xlabel("Read length (nt)")
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/sfigure3/ML_read_length_distribution.png", dpi=300)