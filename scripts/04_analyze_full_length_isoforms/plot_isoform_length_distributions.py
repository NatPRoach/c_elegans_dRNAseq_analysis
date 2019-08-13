#!/usr/bin/env python2

import sys
import sets
import matplotlib
import matplotlib.pyplot as plt

def splitAttributes(attr):
    fields = attr.split(';')
    attr_dict = {}
    for field in fields:
        key,val = field.split('=')
        attr_dict[key] = val
    return attr_dict
def splitGTFAttributes(attr):
    fields = attr.split(';')
    # print attr
    # print fields
    attr_dict = {}
    for field in fields[:-1]:
        key,val = field.strip().split()
        attr_dict[key] = val.strip('\"')
        # print val.strip('\"')
    return attr_dict
infile = open(sys.argv[1])
lengths1 = []
for line in infile:
    fields = line.strip().split()
    length = len(fields[4])
    lengths1.append(length)

# infile = open(sys.argv[2])
# lengths2 = []
# for line in infile:
#     fields = line.strip().split()
#     length = len(fields[4])
#     lengths2.append(length)

splice_isoform_color='#C9B3D5'
annotation_color='#F99B9B'
novel_isoforms_color='#C9B3D5'
font = {"family":"sans-serif",
        "sans-serif":["Helvetica"],
        "size":10}
matplotlib.rc('font',**font)

reference_in = open("../../references/WS265/c_elegans.PRJNA13758.WS265.WormBase.gff3")
isoform_to_length = {}
piRNAs = sets.Set()
ncRNAs = sets.Set()
snoRNAs = sets.Set()
for line in reference_in:
    fields = line.strip().split()
    if fields[1] == "WormBase_transposon":
        continue
    if fields[2] == "exon":
        isoform = fields[8].split(':')[1]
        if isoform in piRNAs or isoform in ncRNAs or isoform in snoRNAs:
            continue
        length = int(fields[4]) - int(fields[3])
        if isoform in isoform_to_length:
            isoform_to_length[isoform] += length
        else:
            isoform_to_length[isoform] = length
    elif fields[2] == "piRNA":
        attr_dict = splitAttributes(fields[8])
        isoform = attr_dict["ID"].split(':')[1]
        piRNAs.add(isoform)
    elif fields[2] == "ncRNA":
        attr_dict = splitAttributes(fields[8])
        isoform = attr_dict["ID"].split(':')[1]
        ncRNAs.add(isoform)
    elif fields[2] == "snoRNA":
        attr_dict = splitAttributes(fields[8])
        isoform = attr_dict["ID"].split(':')[1]
        snoRNAs.add(isoform)
lengths = []
for isoform in isoform_to_length:
    if isoform in piRNAs or isoform in ncRNAs or isoform in snoRNAs:
        continue
    lengths.append(isoform_to_length[isoform])
    # if isoform_to_length[isoform] < 100:
    #     print isoform

# plt.hist(lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.4,label = "annotation isoforms")
# plt.hist(lengths1,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.4,label="our isoforms")

# plt.hist(lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.6,color=annotation_color,edgecolor='black',label = "annotation isoforms",linewidth=0.5)
# plt.hist(lengths1,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.6,color=splice_isoform_color,edgecolor='black',label="our isoforms",linewidth=0.5)
plt.figure(num=None,figsize=(4.0,2.5))

plt.hist(lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.6,edgecolor='black',label = "annotation isoforms",linewidth=0.5)
plt.hist(lengths1,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)

# plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.4,label="our novel isoforms")
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# plt.ylim([0,2500])
plt.ylabel("Fraction of isoforms")
plt.xlabel("Isoform length")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("../../figures/figure1/Figure1E.png",dpi=600)

# plt.show()

# # plt.hist(lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.6,color=annotation_color,edgecolor='black',label = "annotation isoforms",linewidth=0.5)
# # plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.6,color=novel_isoforms_color,edgecolor='black',label="our novel isoforms",linewidth=0.5)
#
# plt.figure(num=None,figsize=(4.0,2.5))
#
# plt.hist(lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.6,edgecolor='black',label = "annotation isoforms",linewidth=0.5)
# plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.6,edgecolor='black',label="all novel isoforms",linewidth=0.5)
# # plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.4,label="our novel isoforms")
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# # plt.ylim([0,2500])
# plt.ylabel("Fraction of isoforms")
# plt.xlabel("Isoform length")
# plt.legend(frameon=False)
# plt.tight_layout()
# plt.savefig("/Users/nproach/Desktop/annotation_vs_our_novel_isoforms.png",dpi=600)
#
# # plt.show()


reference_in = open("../../references/stringtie2/stringtie2_assembly.gtf")
waterston_isoform_to_length = {}
for line in reference_in:
    if line[0] == '#':
        continue
    fields = line.strip().split('\t')
    if fields[2] == "exon":
        attr_dict = splitGTFAttributes(fields[8])
        isoform = attr_dict["transcript_id"]
        length = int(fields[4]) - int(fields[3])
        if isoform in waterston_isoform_to_length:
            waterston_isoform_to_length[isoform] += length
        else:
            waterston_isoform_to_length[isoform] = length
waterston_lengths = []
for isoform in waterston_isoform_to_length:
    # if isoform in piRNAs or isoform in ncRNAs or isoform in snoRNAs:
    #     continue
    waterston_lengths.append(waterston_isoform_to_length[isoform])
    # if isoform_to_length[isoform] < 100:
    #     print isoform

# plt.hist(waterston_lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.6,color=annotation_color,edgecolor='black',label = "waterston isoforms",linewidth=0.5)
# plt.hist(lengths1,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.6,color=splice_isoform_color,edgecolor='black',label="our isoforms",linewidth=0.5)
plt.figure(num=None,figsize=(4.0,2.5))

plt.hist(waterston_lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.6,edgecolor='black',label = "Stringtie2 isoforms",linewidth=0.5)
plt.hist(lengths1,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)

# plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.4,label="our novel isoforms")
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# plt.ylim([0,2500])
plt.ylabel("Fraction of isoforms")
plt.xlabel("Isoform length")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("../../figures/figure1/Figure1F.png",dpi=600)

# # plt.show()
#
# # plt.hist(waterston_lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.6,color=annotation_color,edgecolor='black',label = "waterston isoforms",linewidth=0.5)
# # plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.6,color=novel_isoforms_color,edgecolor='black',label="our novel isoforms",linewidth=0.5)
# plt.figure(num=None,figsize=(4.0,2.5))
#
# plt.hist(waterston_lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.6,edgecolor='black',label = "Stringtie2 isoforms",linewidth=0.5)
# plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.6,edgecolor='black',label="all novel isoforms",linewidth=0.5)
#
# # plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.4,label="our novel isoforms")
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# # plt.ylim([0,2500])
# plt.ylabel("Fraction of isoforms")
# plt.xlabel("Isoform length")
# plt.legend(frameon=False)
# plt.tight_layout()
# plt.savefig("/Users/nproach/Desktop/waterston_vs_novel_isoforms.png",dpi=600)
#
# # plt.show()


# flair_in1 = open("/Users/nproach/Documents/flair/flair.collapse.isoforms.fa")
# flair_lengths = []
# for line in flair_in1:
#     if line[0] == '>':
#         continue
#     flair_lengths.append(len(line.strip()))
#
#
# flair_in2 = open("/Users/nproach/Documents/flair/flair.collapse.isoforms.stringent.fa")
# flair_stringent_lengths = []
# for line in flair_in2:
#     if line[0] == '>':
#         continue
#     flair_stringent_lengths.append(len(line.strip()))


# # plt.hist(waterston_lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.6,color=annotation_color,edgecolor='black',label = "waterston isoforms",linewidth=0.5)
# # plt.hist(lengths1,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.6,color=splice_isoform_color,edgecolor='black',label="our isoforms",linewidth=0.5)
# plt.figure(num=None,figsize=(4.0,2.5))
#
# plt.hist(lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.6,edgecolor='black',label = "annotation isoforms",linewidth=0.5)
# plt.hist(flair_lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.6,edgecolor='black',label="flair isoforms",linewidth=0.5)
#
# # plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.4,label="our novel isoforms")
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# # plt.ylim([0,2500])
# plt.ylabel("Fraction of isoforms")
# plt.xlabel("Isoform length")
# plt.legend()
# plt.tight_layout()
# plt.savefig("/Users/nproach/Desktop/flair_vs_annotation_isoforms.png",dpi=600)
#
# # plt.show()

# # plt.hist(waterston_lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.6,color=annotation_color,edgecolor='black',label = "waterston isoforms",linewidth=0.5)
# # plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.6,color=novel_isoforms_color,edgecolor='black',label="our novel isoforms",linewidth=0.5)
#
# plt.figure(num=None,figsize=(4.0,2.5))
#
# plt.hist(flair_lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.6,edgecolor='black',label = "flair isoforms",linewidth=0.5)
# plt.hist(lengths1,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)
#
# # plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.4,label="our novel isoforms")
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# # plt.ylim([0,2500])
# plt.ylabel("Fraction of isoforms")
# plt.xlabel("Isoform length")
# plt.legend(frameon=False)
# plt.tight_layout()
# plt.savefig("/Users/nproach/Desktop/flair_vs_our_isoforms.png",dpi=600)
#
# # plt.show()

