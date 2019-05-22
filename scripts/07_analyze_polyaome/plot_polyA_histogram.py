#!/usr/bin/env python2
import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np

if len(sys.argv) > 1:
    our_color = sys.argv[1]
    lima_color= sys.argv[2]
else:
    our_color = "#2579B2"
    lima_color = "#389E34"

font = {"family":"sans-serif",
        "sans-serif":["Helvetica"],
        "size":10}
matplotlib.rc('font',**font)
polya_lengths = []
infile2 = open("../../data/L4/bio1/tech1/analysis/ce11_gen_L4_bio1_tech1_Minimap2_k14.polya")
for i,line in enumerate(infile2):
    if i == 0:
        continue
    fields = line.strip().split()
    if fields[9] == "PASS":
        read_len = float(fields[8])
        polya_lengths.append(read_len)
infile2 = open("../../data/L4/bio1/tech2/analysis/ce11_gen_L4_bio1_tech2_Minimap2_k14.polya")
for i,line in enumerate(infile2):
    if i == 0:
        continue
    fields = line.strip().split()
    if fields[9] == "PASS":
        read_len = float(fields[8])
        polya_lengths.append(read_len)
    
infile = open("../../references/polya/polyALengths.txt")
lima_polya_lengths = []
for line in infile:
    lima_polya_lengths.append(float(line.strip()))

plt.figure(num=None,figsize=(3.5,2.5))
sns.kdeplot(lima_polya_lengths,color=lima_color,gridsize=10000,label="Lima et al")
sns.kdeplot(polya_lengths,color = our_color,gridsize=10000,label="This study")
plt.legend(fontsize=10,frameon=False)
plt.xlabel("polyA tail lengths",fontsize=10)
plt.ylabel("Density",fontsize=10)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.xlim(0,200)
plt.tight_layout()
plt.savefig("../../figures/supplementals/sfigure5/sfigure5A.pdf")
