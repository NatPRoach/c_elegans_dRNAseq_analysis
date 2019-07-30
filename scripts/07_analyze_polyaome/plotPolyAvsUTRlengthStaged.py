#!/usr/bin/env python2

import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from scipy import stats

## The following code snippet found at https://zhiyzuo.github.io/Pearson-Correlation-CI-in-Python/
def pearsonr_ci(x,y,alpha=0.05):
    ''' calculate Pearson correlation along with the confidence interval using scipy and numpy
    Parameters
    ----------
    x, y : iterable object such as a list or np.array
      Input for correlation calculation
    alpha : float
      Significance level. 0.05 by default
    Returns
    -------
    r : float
      Pearson's correlation coefficient
    pval : float
      The corresponding p value
    lo, hi : float
      The lower and upper bound of confidence intervals
    '''

    r, p = stats.pearsonr(x,y)
    r_z = np.arctanh(r)
    se = 1/np.sqrt(x.size-3)
    z = stats.norm.ppf(1-alpha/2)
    lo_z, hi_z = r_z-z*se, r_z+z*se
    lo, hi = np.tanh((lo_z, hi_z))
    return r, p, lo, hi

def load_read_id_to_length():
    read_id_to_length = {}### For each read in the polyA file assign its polyA tail length in a dictionary
    tmppolyafile = open("../../data/L1/bio1/tech1/analysis/ce11_gen_L1_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L1/bio1/tech2/analysis/ce11_gen_L1_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L2/bio1/tech1/analysis/ce11_gen_L2_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L2/bio1/tech2/analysis/ce11_gen_L2_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])            

    tmppolyafile = open("../../data/L3/bio1/tech1/analysis/ce11_gen_L3_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L3/bio1/tech2/analysis/ce11_gen_L3_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L4/bio1/tech1/analysis/ce11_gen_L4_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L4/bio1/tech2/analysis/ce11_gen_L4_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])
            
    tmppolyafile = open("../../data/young_adult/bio1/tech1/analysis/ce11_gen_young_adult_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/young_adult/bio1/tech2/analysis/ce11_gen_young_adult_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/adult/bio1/tech1/analysis/ce11_gen_adult_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/adult/bio1/tech2/analysis/ce11_gen_adult_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])
            
    tmppolyafile = open("../../data/male/bio1/tech1/analysis/ce11_gen_male_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/male/bio1/tech2/analysis/ce11_gen_male_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])
    return read_id_to_length


# def plotReadCountsVsPolyA(ax,cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,title):
#     cluster_id_to_read_ids = {}
#     for line in cluster_id_assignments_in:
#         fields = line.strip().split()
#         cluster_id = fields[1]
#         read_id = fields[0]
#         if cluster_id in cluster_id_to_read_ids:
#             cluster_id_to_read_ids[cluster_id].append(read_id)
#         else:
#             cluster_id_to_read_ids[cluster_id] = [read_id]
#
#     median_polya = []
#     log_length = []
#     for cluster_id in cluster_id_to_read_ids:
#         read_ids = cluster_id_to_read_ids[cluster_id]
#         polyas = []
#         counter = 0
#         for read_id in read_ids:
#             if read_id not in read_id_to_length:
#                 continue
#             polyas.append(read_id_to_length[read_id])
#         if len(polyas) >= 10:
#             median_polya.append(np.median(polyas))
#             log_length.append(np.log2(float(cluster_id_to_length[cluster_id])))
#
#     matrix = np.matrix([log_length,median_polya])
#     matrix = matrix.transpose()
#     model = LinearRegression()
#     model = model.fit(matrix[:,0],matrix[:,1])
#     coefs = [model.intercept_[0],model.coef_[0][0]]
#     print coefs
#     x = np.linspace(4,10,1000)
#     y = coefs[0] + coefs[1]*x
#     y_hat = coefs[0] + coefs[1]*np.array(log_length)
#     r2 = r2_score(np.array(median_polya),y_hat)
#     print r2
#     sns.kdeplot(log_length,median_polya,color="#2579B2",ax=ax)
#     ax.plot(x,y,color="#FD7F23")
#     ax.text(4,175,"$R^2$=%.4f"%(r2))
#     ax.set_xlim(3,10)
#     ax.set_ylim(0,200)
#     ax.set_title(title)
#     ax.spines['right'].set_visible(False)
#     ax.spines['top'].set_visible(False)
#     return coefs[1]

def plotReadCountsVsPolyA(ax,cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,title):
    cluster_id_to_read_ids = {}
    for line in cluster_id_assignments_in:
        fields = line.strip().split()
        # if len(fields) != 4:
        #     continue
        cluster_id = fields[1]
        read_id = fields[0]
        #gene_id = fields[1]
        if cluster_id in cluster_id_to_read_ids:
            cluster_id_to_read_ids[cluster_id].append(read_id)
        else:
            cluster_id_to_read_ids[cluster_id] = [read_id]

    median_polya = []
    log_length = []
    # all_polya = []
    # all_expression = []
    for cluster_id in cluster_id_to_read_ids:
        read_ids = cluster_id_to_read_ids[cluster_id]
        polyas = []
        counter = 0
        for read_id in read_ids:
            if read_id not in read_id_to_length:
                continue
            polyas.append(read_id_to_length[read_id])
        if len(polyas) >= 10:
            median_polya.append(np.median(polyas))
            log_length.append(np.log10(float(cluster_id_to_length[cluster_id])))

    matrix = np.matrix([log_length,median_polya])
    matrix = matrix.transpose()
    # poly = PolynomialFeatures(degree=2)
    # model = Pipeline([('poly', PolynomialFeatures(degree=2)),('linear', LinearRegression(fit_intercept=False))])
    model = LinearRegression()
    model = model.fit(matrix[:,0],matrix[:,1])
    coefs = [model.intercept_[0],model.coef_[0][0]]
    print title
    print coefs
    #coefs = model.named_steps['linear'].coef_
    x = np.linspace(1,3,1000)
    y = coefs[0] + coefs[1]*x #+ coefs[:,2]*(x**2)
    y_hat = coefs[0] + coefs[1]*np.array(log_length)
    # r2 = r2_score(np.array(median_polya),y_hat)
    r,p,lo,hi = pearsonr_ci(np.array(median_polya),y_hat)
    r2 = r**2
    print r2
    # print r**2
    print lo**2, hi**2
    print p
    sns.kdeplot(log_length,median_polya,color="#2579B2",ax=ax)
    #ax.scatter(log_length,median_polya,alpha = 0.05,color="#2579B2")
    ax.plot(x,y,color="#FD7F23")
    ax.text(1.05,175,"$R^2$=%.4f"%(r2))
    ax.text(1.05,152.5,"(%.4f,%.4f)"%(lo**2,hi**2))
    ax.text(1.05,120,"p=%e"%(p))
    # ax.ylabel("Median PolyA tail length")
    # ax.xlabel("Log2 Read Counts")
    ax.set_xlim(1,3)
    ax.set_ylim(0,200)
    ax.set_title(title)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    return coefs[1]

def plotSlopes(ax,l1,l2,l3,l4,ya,ga,ml,slope_color="#2579B2"):
    x1 = range(6)
    x2 = range(7)
    y1 = [l1,l2,l3,l4,ya,ga]
    y2 = [l1,l2,l3,l4,ya,ga,ml]
    ax.plot(x1,y1,color=slope_color)
    ax.scatter(x2,y2,color=slope_color)
    ax.set_xticks(range(7))
    ax.set_xticklabels(labels=["L1","L2","L3","L4","yAd","mAd","male"])


font = {"family":"sans-serif",
        "sans-serif":["Helvetica"],
        "size":10}
matplotlib.rc('font',**font)

if len(sys.argv) > 1:
    point_color = sys.argv[1]
    regression_color= sys.argv[2]
else:
    point_color = "#2579B2"
    regression_color = "#FD7F23"

# read_id_to_length = load_read_id_to_length()
# cluster_id_to_length = {}
# bedFile = open("../../results/utrs/beds/all_isoforms_utrs.bed")
# for line in bedFile:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     cluster_id = fields[3]
#     utr_length = end - start
#     cluster_id_to_length[cluster_id] = utr_length
# #plt.figure(figsize=(20,8))
# fig, axes = plt.subplots(2,4,sharex="col",sharey="row",figsize=(8.333,4))
# axes[0,0].set_ylabel("Median PolyA tail length")
# axes[1,0].set_ylabel("Median PolyA tail length")
# axes[1,0].set_xlabel("Log2 3'UTR length")
# axes[1,1].set_xlabel("Log2 3'UTR length")
# axes[1,2].set_xlabel("Log2 3'UTR length")
# axes[1,3].set_xlabel("Log2 3'UTR length")
#
# cluster_id_assignments_in = open("../../results/utrs/assignments/L1_utrs.tsv")
# l1_slope = plotReadCountsVsPolyA(axes[0,0],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"L1")
#
# cluster_id_assignments_in = open("../../results/utrs/assignments/L2_utrs.tsv")
# l2_slope = plotReadCountsVsPolyA(axes[0,1],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"L2")
#
# cluster_id_assignments_in = open("../../results/utrs/assignments/L3_utrs.tsv")
# l3_slope = plotReadCountsVsPolyA(axes[0,2],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"L3")
#
# cluster_id_assignments_in = open("../../results/utrs/assignments/L4_utrs.tsv")
# l4_slope = plotReadCountsVsPolyA(axes[0,3],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"L4")
#
# cluster_id_assignments_in = open("../../results/utrs/assignments/YA_utrs.tsv")
# ya_slope = plotReadCountsVsPolyA(axes[1,0],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"young adult")
#
# cluster_id_assignments_in = open("../../results/utrs/assignments/GA_utrs.tsv")
# ga_slope = plotReadCountsVsPolyA(axes[1,1],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"mature adult")
#
# cluster_id_assignments_in = open("../../results/utrs/assignments/ML_utrs.tsv")
# ml_slope = plotReadCountsVsPolyA(axes[1,2],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"male")
#
# cluster_id_assignments_in = open("../../results/utrs/assignments/all_isoforms_utrs.tsv")
# plotReadCountsVsPolyA(axes[1,3],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"all")
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/sfigure5/sfigure5C.png",dpi=300)
# plt.clf()
#
# fig, axes = plt.subplots(1,2,figsize=(5.0,2),gridspec_kw = {'width_ratios':[1, 2]})
# axes[0].set_ylabel("Median PolyA tail length")
# axes[0].set_xlabel("Log2 3'UTR length")
# axes[1].set_ylabel("Regression line slope")
# axes[0].spines['right'].set_visible(False)
# axes[0].spines['top'].set_visible(False)
# axes[1].spines['right'].set_visible(False)
# axes[1].spines['top'].set_visible(False)
# cluster_id_assignments_in = open("../../results/utrs/assignments/all_isoforms_utrs.tsv")
# plotReadCountsVsPolyA(axes[0],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"")
# plotSlopes(axes[1],l1_slope,l2_slope,l3_slope,l4_slope,ya_slope,ga_slope,ml_slope,regression_color)
# plt.tight_layout()
# plt.savefig("../../figures/figure4/figure4D_left_and_center.png",dpi=450)

### Stringent


read_id_to_length = load_read_id_to_length()
cluster_id_to_length = {}
bedFile = open("../../results/utrs/beds/all_stringent_isoforms_utrs.bed")
for line in bedFile:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    cluster_id = fields[3]
    utr_length = end - start
    cluster_id_to_length[cluster_id] = utr_length
#plt.figure(figsize=(20,8))
fig, axes = plt.subplots(2,4,sharex="col",sharey="row",figsize=(8.333,4))
axes[0,0].set_ylabel("Median PolyA tail length")
axes[1,0].set_ylabel("Median PolyA tail length")
axes[1,0].set_xlabel("Log10 3'UTR length")
axes[1,1].set_xlabel("Log10 3'UTR length")
axes[1,2].set_xlabel("Log10 3'UTR length")
axes[1,3].set_xlabel("Log10 3'UTR length")

cluster_id_assignments_in = open("../../results/utrs/assignments/L1_stringent_utrs.tsv")
l1_slope = plotReadCountsVsPolyA(axes[0,0],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"L1")

cluster_id_assignments_in = open("../../results/utrs/assignments/L2_stringent_utrs.tsv")
l2_slope = plotReadCountsVsPolyA(axes[0,1],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"L2")

cluster_id_assignments_in = open("../../results/utrs/assignments/L3_stringent_utrs.tsv")
l3_slope = plotReadCountsVsPolyA(axes[0,2],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"L3")

cluster_id_assignments_in = open("../../results/utrs/assignments/L4_stringent_utrs.tsv")
l4_slope = plotReadCountsVsPolyA(axes[0,3],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"L4")

cluster_id_assignments_in = open("../../results/utrs/assignments/YA_stringent_utrs.tsv")
ya_slope = plotReadCountsVsPolyA(axes[1,0],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"young adult")

cluster_id_assignments_in = open("../../results/utrs/assignments/GA_stringent_utrs.tsv")
ga_slope = plotReadCountsVsPolyA(axes[1,1],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"mature adult")

cluster_id_assignments_in = open("../../results/utrs/assignments/ML_stringent_utrs.tsv")
ml_slope = plotReadCountsVsPolyA(axes[1,2],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"male")

cluster_id_assignments_in = open("../../results/utrs/assignments/all_stringent_isoforms_utrs.tsv")
plotReadCountsVsPolyA(axes[1,3],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"all")
plt.tight_layout()
plt.savefig("../../figures/supplementals/sfigure5/sfigure5C_stringent.png",dpi=300)
plt.clf()

fig, axes = plt.subplots(1,2,figsize=(5.0,2),gridspec_kw = {'width_ratios':[1, 2]})
axes[0].set_ylabel("Median PolyA tail length")
axes[0].set_xlabel("Log10 3'UTR length")
axes[1].set_ylabel("Regression line slope")
axes[0].spines['right'].set_visible(False)
axes[0].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[1].spines['top'].set_visible(False)
cluster_id_assignments_in = open("../../results/utrs/assignments/all_stringent_isoforms_utrs.tsv")
plotReadCountsVsPolyA(axes[0],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"")
plotSlopes(axes[1],l1_slope,l2_slope,l3_slope,l4_slope,ya_slope,ga_slope,ml_slope,regression_color)
plt.tight_layout()
plt.savefig("../../figures/figure4/figure4D_left_and_center_stringent.png",dpi=450)


### Sensitive

read_id_to_length = load_read_id_to_length()
cluster_id_to_length = {}
bedFile = open("../../results/utrs/beds/all_sensitive_isoforms_utrs.bed")
for line in bedFile:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    cluster_id = fields[3]
    utr_length = end - start
    cluster_id_to_length[cluster_id] = utr_length
#plt.figure(figsize=(20,8))
fig, axes = plt.subplots(2,4,sharex="col",sharey="row",figsize=(8.333,4))
axes[0,0].set_ylabel("Median PolyA tail length")
axes[1,0].set_ylabel("Median PolyA tail length")
axes[1,0].set_xlabel("Log10 3'UTR length")
axes[1,1].set_xlabel("Log10 3'UTR length")
axes[1,2].set_xlabel("Log10 3'UTR length")
axes[1,3].set_xlabel("Log10 3'UTR length")

cluster_id_assignments_in = open("../../results/utrs/assignments/L1_sensitive_utrs.tsv")
l1_slope = plotReadCountsVsPolyA(axes[0,0],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"L1")

cluster_id_assignments_in = open("../../results/utrs/assignments/L2_sensitive_utrs.tsv")
l2_slope = plotReadCountsVsPolyA(axes[0,1],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"L2")

cluster_id_assignments_in = open("../../results/utrs/assignments/L3_sensitive_utrs.tsv")
l3_slope = plotReadCountsVsPolyA(axes[0,2],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"L3")

cluster_id_assignments_in = open("../../results/utrs/assignments/L4_sensitive_utrs.tsv")
l4_slope = plotReadCountsVsPolyA(axes[0,3],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"L4")

cluster_id_assignments_in = open("../../results/utrs/assignments/YA_sensitive_utrs.tsv")
ya_slope = plotReadCountsVsPolyA(axes[1,0],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"young adult")

cluster_id_assignments_in = open("../../results/utrs/assignments/GA_sensitive_utrs.tsv")
ga_slope = plotReadCountsVsPolyA(axes[1,1],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"mature adult")

cluster_id_assignments_in = open("../../results/utrs/assignments/ML_sensitive_utrs.tsv")
ml_slope = plotReadCountsVsPolyA(axes[1,2],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"male")

cluster_id_assignments_in = open("../../results/utrs/assignments/all_sensitive_isoforms_utrs.tsv")
plotReadCountsVsPolyA(axes[1,3],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"all")
plt.tight_layout()
plt.savefig("../../figures/supplementals/sfigure5/sfigure5C_sensitive.png",dpi=300)
plt.clf()

fig, axes = plt.subplots(1,2,figsize=(5.0,2),gridspec_kw = {'width_ratios':[1, 2]})
axes[0].set_ylabel("Median PolyA tail length")
axes[0].set_xlabel("Log10 3'UTR length")
axes[1].set_ylabel("Regression line slope")
axes[0].spines['right'].set_visible(False)
axes[0].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[1].spines['top'].set_visible(False)
cluster_id_assignments_in = open("../../results/utrs/assignments/all_sensitive_isoforms_utrs.tsv")
plotReadCountsVsPolyA(axes[0],cluster_id_assignments_in,read_id_to_length,cluster_id_to_length,"")
plotSlopes(axes[1],l1_slope,l2_slope,l3_slope,l4_slope,ya_slope,ga_slope,ml_slope,regression_color)
plt.tight_layout()
plt.savefig("../../figures/figure4/figure4D_left_and_center_sensitive.png",dpi=450)