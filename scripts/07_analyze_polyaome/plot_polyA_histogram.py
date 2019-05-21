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
# infile = open("limaEtAlPolyALength.txt")
# plt.set_cmap("Paired")
# x = []
# y = []
# xscale_ratios = []
# yscale_ratios = []
# # xscale_tuples = []
# # yscale_tuples = []
# for line in infile:
#     fields = line.strip().split()
#     x.append(float(fields[0]))
#     y.append(float(fields[1]))
#     if fields[2] != '-':
#         xscale_ratios.append(float(fields[2])/float(fields[0]))
#     if fields[3] != '-':
#         yscale_ratios.append(float(fields[3])/float(fields[1]))
#
# xscale = np.mean(xscale_ratios)
# # yscale = np.mean(yscale_ratios)
# x_prime = []
# for i in range(len(x)):
#     x_prime.append(x[i]*xscale)
#
# auc = 0.
# for i in range(len(x)-1):
#     height = (y[i+1] + y[i]) / 2
#     base = x_prime[i+1] - x_prime[i]
#     auc += height * base
# y_prime = []
# for i in range(len(x)):
#     y_prime.append(y[i]/auc)
#
# auc3 = 0.
# for i in range(len(x)-1):
#     height = (y_prime[i+1] + y_prime[i]) / 2
#     base = x_prime[i+1] - x_prime[i]
#     auc3 += height * base
# print auc3
    # y_prime.append(y[i]*yscale)
# density_test = []
# for i in range(len(x_prime)):
#     x_i = x_prime[i]
#     for j in range(int(y[i])):
#         density_test.append(x_i)
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
    
# sns.distplot(density_test)
# sns.distplot(polya_lengths,bins=range(1000))
# ax = sns.kdeplot(polya_lengths)
# x_vals,y_vals = ax.lines[0].get_data()
# plt.clf()
# auc2 = 0.
# for i in range(len(x_vals)-1):
#     height = (y_vals[i+1] + y_vals[i]) / 2
#     base = x_vals[i+1] - x_vals[i]
#     auc2 += height * base
# print auc2
#
# y_vals2 = []
# for i in range(len(x_vals)):
#     y_vals2.append(y_vals[i]/auc2)
infile = open("../../references/polya/polyALengths.txt")
lima_polya_lengths = []
for line in infile:
    lima_polya_lengths.append(float(line.strip()))

# plt.hist(polya_lengths,density=True,bins = range(200),color="#2579B2",alpha=0.25,label="Our L4 poly(A) lengths")
# plt.hist(lima_polya_lengths,density=True,bins = range(200),color="#389E34",alpha=0.25,label="Lima et al poly(A) lengths")

# plt.plot(x_vals,y_vals2)
#plt.plot(x_prime,y_prime,color="#FD7F23",linewidth=2,label="Lima et al poly(A) lengths")
plt.figure(num=None,figsize=(3.5,2.5))
sns.kdeplot(lima_polya_lengths,color=lima_color,gridsize=10000,label="Lima et al")
sns.kdeplot(polya_lengths,color = our_color,gridsize=10000,label="This study")
plt.legend(fontsize=10,frameon=False)
plt.xlabel("polyA tail lengths",fontsize=10)
plt.ylabel("Density",fontsize=10)
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
# n, bins, rectangles = plt.hist(polya_lengths,density=True,bins = range(200))
# # plt.plot(x_vals,y_vals2)
# plt.plot(x_prime,y_prime)
plt.xlim(0,200)
plt.tight_layout()
plt.savefig("../../figures/supplementals/sfigure5/sfigure5A.pdf")
# plt.show()