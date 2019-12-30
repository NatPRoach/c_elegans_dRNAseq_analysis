#!/usr/bin/env python2

import sys
import sets
import matplotlib
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
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

def getFCInWindows(lengths1,lengths2=[],window_size=1000,start=0,end=6000):
    bins = range(start,end + 1,window_size)
    counts1 = [0] * (len(bins) - 1)
    counts2 = [0] * (len(bins) - 1)
    lengths1_size = len(lengths1)
    lengths2_size = len(lengths2)
    for length in lengths1:
        lo = 0
        hi = len(bins)
        i = len(bins) / 2
        flag = bins[i] <= length < bins[i+1]
        while not flag:
            if bins[i] > length:
                hi = i
            elif bins[i+1] <= length:
                lo = i+1
            i = (hi + lo) / 2
            if i < len(bins) - 1:
                flag = bins[i] <= length < bins[i+1]
            else:
                i = -1
                break
        if i != -1:
            counts1[i] += 1
    for length in lengths2:
        lo = 0
        hi = len(bins)
        i = len(bins) / 2
        flag = bins[i] <= length < bins[i+1]
        while not flag:
            if bins[i] > length:
                hi = i
            elif bins[i+1] <= length:
                lo = i+1
            i = (hi + lo) / 2
            if i < len(bins) - 1:
                flag = bins[i] <= length < bins[i+1]
            else:
                i = -1
                break
        if i != -1:
            counts2[i] += 1
    fcs = []
    for i in range(len(counts1)):
        # fcs.append(np.log2((float(counts1[i]) / float(lengths1_size))/ (float(counts2[i])/ float(lengths2_size ))))
        fcs.append(((float(counts1[i]) / float(lengths1_size))/ (float(counts2[i])/ float(lengths2_size ))))
    # print fcs
    return fcs
    

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
# plt.figure(num=None,figsize=(4.0,2.5))
fig,[ax1,ax2] = plt.subplots(2,figsize=(4.0,2.5),sharex=True,gridspec_kw = {'height_ratios':[2, 1]})

# plt.hist(lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.6,edgecolor='black',label = "annotation isoforms",linewidth=0.5)
# plt.hist(lengths1,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)
fc = getFCInWindows(lengths1,lengths)
ax2.plot([0,6000],[1,1],color="red",linestyle="--")
for i in range(len(fc)):
    # x_location = i * 1000 + 500
    x1 = i * 1000
    x2 = ((i+1) * 1000) - 1
    x = [x1,x2]
    y = [fc[i],fc[i]]
    ax2.add_patch(mpatches.Rectangle((x1,0),width = 1000,height= fc[i],zorder = 0,edgecolor="black"))

# y_location = 0.0007
# plt.text(-650,y_location,"FC:",horizontalalignment='center')
# for i in range(len(fc)):
#     x_location = i * 1000 + 500
#     plt.text(x_location,y_location,"%0.3f" %(fc[i]),horizontalalignment='center')
#     if i != len(fc) - 1:
#         x_location = (i + 1) *  1000
#         plt.axvline(x_location,ymin=0.9,c='gray',linestyle='-')


# plt.hist(lengths,bins = [x - 0.5 for x in range(0,40000,100)],density=True,alpha = 0.6,edgecolor='black',label = "annotation isoforms",linewidth=0.5)
# plt.hist(lengths1,bins = [x - 0.5 for x in range(0,40000,100)],density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)
ax1.hist(lengths,bins = [x - 0.5 for x in range(0,40000,100)],density=True,alpha = 0.6,edgecolor='black',label = "annotation isoforms",linewidth=0.5)
ax1.hist(lengths1,bins = [x - 0.5 for x in range(0,40000,100)],density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)


# plt.hist(lengths,bins = 400,density=True,alpha = 0.6,edgecolor='black',label = "annotation isoforms",linewidth=0.5)
# plt.hist(lengths1,bins = 400,density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)
#

# plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.4,label="our novel isoforms")
ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
# plt.ylim([0,2500])
ax1.set_ylim(0,0.00075)
ax1.set_ylabel("Fraction of\nisoforms")
ax1.set_yticks([0,0.00025,0.0005,0.00075])
ax1.set_yticklabels(["0","2.5e-4","5e-4","7.5e-4"])
# ax2.set_ylim(-0.5, 1.25)
ax2.set_ylim(0.8, 1.25)
ax2.set_ylabel("dRNA / Annotation\nFold Change")
ax2.set_yticks([0.8,1,1.2])
plt.xlim([0,6000])
plt.xlabel("Isoform length")
ax1.legend(frameon=False)
# ax = plt.gca()

# # plt.ylim([0,2500])
# plt.xlim([0,6000])
# plt.ylim(0,0.00075)
# plt.ylabel("Fraction of isoforms")
# plt.xlabel("Isoform length")
# plt.legend(facecolor="white",framealpha=1.0,edgecolor="white",loc=[0.3,0.5])
plt.tight_layout()
plt.savefig("../../figures/figure1/figure1E.png",dpi=600)

# plt.figure(num=None,figsize=(4.0,2.5))
#
# # plt.hist(lengths,bins = [x - 0.5 for x in range(0,6000,100)],density=True,alpha = 0.6,edgecolor='black',label = "annotation isoforms",linewidth=0.5)
# # plt.hist(lengths1,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)
#
# plt.hist(lengths,bins = [x - 0.5 for x in range(0,40000,100)],density=True,alpha = 0.6,edgecolor='black',label = "annotation isoforms",linewidth=0.5)
# plt.hist(lengths1,bins = [x - 0.5 for x in range(0,40000,100)],density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)
#
#
# # plt.hist(lengths,bins = 400,density=True,alpha = 0.6,edgecolor='black',label = "annotation isoforms",linewidth=0.5)
# # plt.hist(lengths1,bins = 400,density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)
# #
#
# # plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.4,label="our novel isoforms")
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# # plt.ylim([0,2500])
# plt.xlim([6000,20000])
# plt.ylim([0,0.000015])
# plt.ylabel("Fraction of isoforms")
# plt.xlabel("Isoform length")
# plt.legend(frameon=False)
# plt.tight_layout()
# plt.savefig("../../figures/figure1/figure1E_6000-10000.png",dpi=600)
#
print np.mean(lengths)
print np.mean(lengths1)

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
# plt.figure(num=None,figsize=(4.0,2.5))
fig,[ax1,ax2] = plt.subplots(2,figsize=(4.0,2.5),sharex=True,gridspec_kw = {'height_ratios':[2, 1]})
# plt.hist(waterston_lengths,bins = 400,density=True,alpha = 0.6,edgecolor='black',label = "Stringtie2 isoforms",linewidth=0.5)
# plt.hist(lengths1,bins = 400,density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)

# ax1 = plt.gca()
fc = getFCInWindows(lengths1,waterston_lengths)
# ax2 = ax1.twinx()

ax2.plot([0,6000],[1,1],color="red",linestyle="--")
# y_location = 0.0007
# plt.text(-650,y_location,"FC:",horizontalalignment='center')
for i in range(len(fc)):
    # x_location = i * 1000 + 500
    x1 = i * 1000
    x2 = ((i+1) * 1000) - 1
    x = [x1,x2]
    y = [fc[i],fc[i]]
    ax2.add_patch(mpatches.Rectangle((x1,0),width = 1000,height= fc[i],zorder = 0,edgecolor="black"))
    # ax2.plot(x,y,color="red",zorder = (i+1) * 5)
    # plt.text(x_location,y_location,"%0.3f" %(fc[i]),horizontalalignment='center')
    # if i != len(fc) - 1:
    #     x_location = (i + 1) *  1000
    #     plt.axvline(x_location,ymin=0.9,c='gray',linestyle='-')
# plt.hist(waterston_lengths,bins = [x - 0.5 for x in range(0,40000,100)],density=True,alpha = 0.6,edgecolor='black',label = "Stringtie2 isoforms",linewidth=0.5)
# plt.hist(lengths1,bins = [x - 0.5 for x in range(0,40000,100)],density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)
# plt.plot(-1,-1,color="red",label = "Fold Changes")
ax1.hist(waterston_lengths,bins = [x - 0.5 for x in range(0,40000,100)],density=True,alpha = 0.6,edgecolor='black',label = "StringTie2 isoforms",linewidth=0.5,zorder=100)
ax1.hist(lengths1,bins = [x - 0.5 for x in range(0,40000,100)],density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5,zorder = 105)
# ax1.plot(-1,-1,color="red",label = "Fold Changes")


# plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.4,label="our novel isoforms")

ax1.spines['right'].set_visible(False)
ax1.spines['top'].set_visible(False)
# ax2.spines['right'].set_visible(False)
# ax2.spines['top'].set_visible(False)




# plt.ylim([0,2500])
ax1.set_ylim(0,0.00075)
ax1.set_ylabel("Fraction of\nisoforms")
ax1.set_yticks([0,0.00025,0.0005,0.00075])
ax1.set_yticklabels(["0","2.5e-4","5e-4","7.5e-4"])
# ax2.set_ylim(-0.5, 1.25)
ax2.set_ylim(0.6, 1.25)
ax2.set_ylabel("dRNA / StringTie2\nFold Change")
ax2.set_yticks([0.6,0.8,1,1.2])
plt.xlim([0,6000])
plt.xlabel("Isoform length")
ax1.legend(frameon=False)
# ax1.legend()
# ax.legend(loc="center right",frameon=False)
plt.tight_layout()
plt.savefig("../../figures/figure1/figure1F.png",dpi=600)


plt.figure(num=None,figsize=(4.0,2.5))

# plt.hist(waterston_lengths,bins = 400,density=True,alpha = 0.6,edgecolor='black',label = "Stringtie2 isoforms",linewidth=0.5)
# plt.hist(lengths1,bins = 400,density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)

plt.hist(lengths,bins = [x - 0.5 for x in range(0,40000,100)],density=True,alpha = 0.6,edgecolor='black',label = "annotation isoforms",linewidth=0.5)
plt.hist(lengths1,bins = [x - 0.5 for x in range(0,40000,100)],density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)

# plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.4,label="our novel isoforms")
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# plt.ylim([0,2500])
plt.xlim([6000,20000])
plt.ylim([0,0.000015])
plt.ylabel("Fraction of isoforms")
plt.xlabel("Isoform length")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("../../figures/supplementals/sfigure2/sfigure2C.png",dpi=600)

plt.clf()

plt.figure(num=None,figsize=(4.0,2.5))

# plt.hist(waterston_lengths,bins = 400,density=True,alpha = 0.6,edgecolor='black',label = "Stringtie2 isoforms",linewidth=0.5)
# plt.hist(lengths1,bins = 400,density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)

plt.hist(waterston_lengths,bins = [x - 0.5 for x in range(0,40000,100)],density=True,alpha = 0.6,edgecolor='black',label = "StringTie2 isoforms",linewidth=0.5)
plt.hist(lengths1,bins = [x - 0.5 for x in range(0,40000,100)],density=True, alpha = 0.6,edgecolor='black',label="all isoforms",linewidth=0.5)

# plt.hist(lengths2,bins = [x - 0.5 for x in range(0,6000,100)],density=True, alpha = 0.4,label="our novel isoforms")
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# plt.ylim([0,2500])
plt.xlim([6000,20000])
plt.ylim([0,0.000015])
plt.ylabel("Fraction of isoforms")
plt.xlabel("Isoform length")
plt.legend(frameon=False)
plt.tight_layout()
plt.savefig("../../figures/supplementals/sfigure2/sfigure2D.png",dpi=600)

plt.clf()

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

