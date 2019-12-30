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
from statsmodels.stats.multitest import multipletests 
#from sklearn.pipeline import Pipeline

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

def plotReadCountsVsPolyA(ax,tx_id_assignments_in,read_id_to_length,title,point_color="#2579B2",regression_color="#FD7F23"):
    tx_id_to_read_ids = {}
    for line in tx_id_assignments_in:
        fields = line.strip().split()
        if len(fields) != 4:
            continue
        tx_id = fields[3]
        read_id = fields[0]
        gene_id = fields[1]
        if tx_id in tx_id_to_read_ids:
            tx_id_to_read_ids[tx_id].append(read_id)
        else:
            tx_id_to_read_ids[tx_id] = [read_id]

    median_polya = []
    log_expression = []
    all_polya = []
    all_expression = []
    for tx_id in tx_id_to_read_ids:
        read_ids = tx_id_to_read_ids[tx_id]
        polyas = []
        counter = 0
        for read_id in read_ids:
            if read_id not in read_id_to_length:
                continue
            polyas.append(read_id_to_length[read_id])
        if len(polyas) > 9:
            median_polya.append(np.median(polyas))
            log_expression.append(np.log10(float(len(read_ids))))

    matrix = np.matrix([log_expression,median_polya])
    matrix = matrix.transpose()
    # poly = PolynomialFeatures(degree=2)
    # model = Pipeline([('poly', PolynomialFeatures(degree=2)),('linear', LinearRegression(fit_intercept=False))])
    model = LinearRegression()
    model = model.fit(matrix[:,0],matrix[:,1])
    coefs = [model.intercept_[0],model.coef_[0][0]]
    print title
    print coefs
    #coefs = model.named_steps['linear'].coef_
    # x = np.linspace(2.3,12,1000)
    x = np.linspace(0,5,1000)
    y = coefs[0] + coefs[1]*x #+ coefs[:,2]*(x**2)
    y_hat = coefs[0] + coefs[1]*np.array(log_expression)
    # r2 = r2_score(np.array(median_polya),y_hat)
    # print r2
    r,p,lo,hi = pearsonr_ci(np.array(median_polya),y_hat)
    r2 = r**2
    print r2
    # print r**2
    print lo**2, hi**2
    print p
    
    # ax.scatter(log_expression,median_polya,alpha = 0.05,color=point_color)
    sns.kdeplot(log_expression,median_polya,color=point_color,ax=ax, n_levels=15)
    ax.plot(x,y,color=regression_color)
    # ax.text(2.3,175,"$R^2$=%.4f"%(r2))
    ax.text(2.05 -.5,175+10,"$R^2$=%.4f"%(r2))
    ax.text(2.05-.5,152.5+10,"(%.4f,%.4f)"%(lo**2,hi**2))
    # ax.text(2.05-.5,120+10,"p=%.4e"%(p))
    # ax.ylabel("Median PolyA tail length")
    # ax.xlabel("Log2 Read Counts")
    ax.set_xlim((0,5))
    # ax.set_xlim((1.3,16))
    ax.set_ylim(0,200)
    ax.set_title(title)
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    return coefs[1],p

def plotPval(ax,p):
    ax.text(2.05-.5,120+10,"p=%.4e"%(p))

def plotReadCountsVsPolyALeftPanel(ax,tx_id_assignments_in,read_id_to_length,title,point_color="#2579B2",regression_color="#FD7F23"):
    tx_id_to_read_ids = {}
    for line in tx_id_assignments_in:
        fields = line.strip().split()
        if len(fields) != 4:
            continue
        tx_id = fields[3]
        read_id = fields[0]
        gene_id = fields[1]
        if tx_id in tx_id_to_read_ids:
            tx_id_to_read_ids[tx_id].append(read_id)
        else:
            tx_id_to_read_ids[tx_id] = [read_id]

    median_polya = []
    log_expression = []
    all_polya = []
    all_expression = []
    for tx_id in tx_id_to_read_ids:
        read_ids = tx_id_to_read_ids[tx_id]
        polyas = []
        counter = 0
        for read_id in read_ids:
            if read_id not in read_id_to_length:
                continue
            polyas.append(read_id_to_length[read_id])
        if len(polyas) > 9:
            median_polya.append(np.median(polyas))
            log_expression.append(np.log10(float(len(read_ids))))

    matrix = np.matrix([log_expression,median_polya])
    matrix = matrix.transpose()
    # poly = PolynomialFeatures(degree=2)
    # model = Pipeline([('poly', PolynomialFeatures(degree=2)),('linear', LinearRegression(fit_intercept=False))])
    model = LinearRegression()
    model = model.fit(matrix[:,0],matrix[:,1])
    coefs = [model.intercept_[0],model.coef_[0][0]]
    print title
    print coefs
    #coefs = model.named_steps['linear'].coef_
    # x = np.linspace(2.3,12,1000)
    x = np.linspace(0,5,1000)
    y = coefs[0] + coefs[1]*x #+ coefs[:,2]*(x**2)
    y_hat = coefs[0] + coefs[1]*np.array(log_expression)
    # r2 = r2_score(np.array(median_polya),y_hat)
    # print r2
    r,p,lo,hi = pearsonr_ci(np.array(median_polya),y_hat)
    r2 = r**2
    print r2
    # print r**2
    print lo**2, hi**2
    print p
    
    # ax.scatter(log_expression,median_polya,alpha = 0.05,color=point_color)
    sns.kdeplot(log_expression,median_polya,color=point_color,ax=ax, n_levels=15)
    ax.plot(x,y,color=regression_color)
    # ax.text(2.3,175,"$R^2$=%.4f"%(r2))
    ax.text(2.05 -.5,175+10,"$R^2$=%.4f"%(r2))
    # ax.text(2.05-.5,152.5+10,"(%.4f,%.4f)"%(lo**2,hi**2))
    # ax.text(2.05-.5,120+10,"p=%.4e"%(p))
    # ax.ylabel("Median PolyA tail length")
    # ax.xlabel("Log2 Read Counts")
    ax.set_xlim((0,5))
    # ax.set_xlim((1.3,16))
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
# fig, axes = plt.subplots(2,4,sharex="col",sharey="row",figsize=(8.333,4))
# axes[0,0].set_ylabel("Median PolyA tail length")
# axes[1,0].set_ylabel("Median PolyA tail length")
# axes[1,0].set_xlabel("Log Read Counts")
# axes[1,1].set_xlabel("Log Read Counts")
# axes[1,2].set_xlabel("Log Read Counts")
# axes[1,3].set_xlabel("Log Read Counts")
#
# tx_id_assignments_in = open("../../results/correctionLogs/L1.gene.txt")
# l1_slope = plotReadCountsVsPolyA(axes[0,0],tx_id_assignments_in,read_id_to_length,"L1")
#
# tx_id_assignments_in = open("../../results/correctionLogs/L2.gene.txt")
# l2_slope = plotReadCountsVsPolyA(axes[0,1],tx_id_assignments_in,read_id_to_length,"L2")
#
# tx_id_assignments_in = open("../../results/correctionLogs/L3.gene.txt")
# l3_slope = plotReadCountsVsPolyA(axes[0,2],tx_id_assignments_in,read_id_to_length,"L3")
#
# tx_id_assignments_in = open("../../results/correctionLogs/L4.gene.txt")
# l4_slope = plotReadCountsVsPolyA(axes[0,3],tx_id_assignments_in,read_id_to_length,"L4")
#
# tx_id_assignments_in = open("../../results/correctionLogs/YA.gene.txt")
# ya_slope = plotReadCountsVsPolyA(axes[1,0],tx_id_assignments_in,read_id_to_length,"young adult")
#
# tx_id_assignments_in = open("../../results/correctionLogs/GA.gene.txt")
# ga_slope = plotReadCountsVsPolyA(axes[1,1],tx_id_assignments_in,read_id_to_length,"mature adult")
#
# tx_id_assignments_in = open("../../results/correctionLogs/ML.gene.txt")
# ml_slope = plotReadCountsVsPolyA(axes[1,2],tx_id_assignments_in,read_id_to_length,"male")
# tx_id_assignments_in = open("../../results/correctionLogs/all.gene.txt")
# plotReadCountsVsPolyA(axes[1,3],tx_id_assignments_in,read_id_to_length,"all")
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/sfigure5/sfigure5B.png",dpi=450)
# plt.clf()
#
# fig, axes = plt.subplots(1,2,figsize=(5.0,2),gridspec_kw = {'width_ratios':[1, 2]})
# axes[0].set_ylabel("Median PolyA tail length")
# axes[0].set_xlabel("Log Read Counts")
# axes[1].set_ylabel("Regression line slope")
# axes[0].spines['right'].set_visible(False)
# axes[0].spines['top'].set_visible(False)
# axes[1].spines['right'].set_visible(False)
# axes[1].spines['top'].set_visible(False)
# tx_id_assignments_in = open("../../results/correctionLogs/all.gene.txt")
# plotReadCountsVsPolyA(axes[0],tx_id_assignments_in,read_id_to_length,"")
# plotSlopes(axes[1],l1_slope,l2_slope,l3_slope,l4_slope,ya_slope,ga_slope,ml_slope,regression_color)
# plt.tight_layout()
# plt.savefig("../../figures/figure4/figure4C_left_and_center.png",dpi=450)


# ### Sensitive
#
# read_id_to_length = load_read_id_to_length()
# fig, axes = plt.subplots(2,4,sharex="col",sharey="row",figsize=(8.333,4))
# axes[0,0].set_ylabel("Median PolyA tail length")
# axes[1,0].set_ylabel("Median PolyA tail length")
# axes[1,0].set_xlabel("Log10 Read Counts")
# axes[1,1].set_xlabel("Log10 Read Counts")
# axes[1,2].set_xlabel("Log10 Read Counts")
# axes[1,3].set_xlabel("Log10 Read Counts")
#
# tx_id_assignments_in = open("../../results/correctionLogs/L1_sensitive.gene.txt")
# l1_slope = plotReadCountsVsPolyA(axes[0,0],tx_id_assignments_in,read_id_to_length,"L1")
#
# tx_id_assignments_in = open("../../results/correctionLogs/L2_sensitive.gene.txt")
# l2_slope = plotReadCountsVsPolyA(axes[0,1],tx_id_assignments_in,read_id_to_length,"L2")
#
# tx_id_assignments_in = open("../../results/correctionLogs/L3_sensitive.gene.txt")
# l3_slope = plotReadCountsVsPolyA(axes[0,2],tx_id_assignments_in,read_id_to_length,"L3")
#
# tx_id_assignments_in = open("../../results/correctionLogs/L4_sensitive.gene.txt")
# l4_slope = plotReadCountsVsPolyA(axes[0,3],tx_id_assignments_in,read_id_to_length,"L4")
#
# tx_id_assignments_in = open("../../results/correctionLogs/YA_sensitive.gene.txt")
# ya_slope = plotReadCountsVsPolyA(axes[1,0],tx_id_assignments_in,read_id_to_length,"young adult")
#
# tx_id_assignments_in = open("../../results/correctionLogs/GA_sensitive.gene.txt")
# ga_slope = plotReadCountsVsPolyA(axes[1,1],tx_id_assignments_in,read_id_to_length,"mature adult")
#
# tx_id_assignments_in = open("../../results/correctionLogs/ML_sensitive.gene.txt")
# ml_slope = plotReadCountsVsPolyA(axes[1,2],tx_id_assignments_in,read_id_to_length,"male")
# tx_id_assignments_in = open("../../results/correctionLogs/all_sensitive.gene.txt")
# plotReadCountsVsPolyA(axes[1,3],tx_id_assignments_in,read_id_to_length,"all")
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/sfigure5/sfigure5B_sensitive.png",dpi=450)
# plt.clf()
#
# fig, axes = plt.subplots(1,2,figsize=(5.0,2),gridspec_kw = {'width_ratios':[1, 2]})
# axes[0].set_ylabel("Median PolyA tail length")
# axes[0].set_xlabel("Log10 Read Counts")
# axes[1].set_ylabel("Regression line slope")
# axes[0].spines['right'].set_visible(False)
# axes[0].spines['top'].set_visible(False)
# axes[1].spines['right'].set_visible(False)
# axes[1].spines['top'].set_visible(False)
# tx_id_assignments_in = open("../../results/correctionLogs/all_sensitive.gene.txt")
# plotReadCountsVsPolyA(axes[0],tx_id_assignments_in,read_id_to_length,"")
# plotSlopes(axes[1],l1_slope,l2_slope,l3_slope,l4_slope,ya_slope,ga_slope,ml_slope,regression_color)
# plt.tight_layout()
# plt.savefig("../../figures/figure4/figure4C_left_and_center_sensitive.png",dpi=450)




### Stringent

read_id_to_length = load_read_id_to_length()
fig, axes = plt.subplots(2,4,sharex="col",sharey="row",figsize=(8.333,4))
axes[0,0].set_ylabel("Median Poly(A) tail length")
axes[1,0].set_ylabel("Median Poly(A) tail length")
axes[1,0].set_xlabel("Log$_{10}$ Read Counts")
axes[1,1].set_xlabel("Log$_{10}$ Read Counts")
axes[1,2].set_xlabel("Log$_{10}$ Read Counts")
axes[1,3].set_xlabel("Log$_{10}$ Read Counts")

tx_id_assignments_in = open("../../results/correctionLogs/L1_stringent.gene.txt")
l1_slope,l1_pval = plotReadCountsVsPolyA(axes[0,0],tx_id_assignments_in,read_id_to_length,"L1")

tx_id_assignments_in = open("../../results/correctionLogs/L2_stringent.gene.txt")
l2_slope,l2_pval = plotReadCountsVsPolyA(axes[0,1],tx_id_assignments_in,read_id_to_length,"L2")

tx_id_assignments_in = open("../../results/correctionLogs/L3_stringent.gene.txt")
l3_slope,l3_pval = plotReadCountsVsPolyA(axes[0,2],tx_id_assignments_in,read_id_to_length,"L3")

tx_id_assignments_in = open("../../results/correctionLogs/L4_stringent.gene.txt")
l4_slope,l4_pval = plotReadCountsVsPolyA(axes[0,3],tx_id_assignments_in,read_id_to_length,"L4")

tx_id_assignments_in = open("../../results/correctionLogs/YA_stringent.gene.txt")
ya_slope,ya_pval = plotReadCountsVsPolyA(axes[1,0],tx_id_assignments_in,read_id_to_length,"young adult")

tx_id_assignments_in = open("../../results/correctionLogs/GA_stringent.gene.txt")
ga_slope,ga_pval = plotReadCountsVsPolyA(axes[1,1],tx_id_assignments_in,read_id_to_length,"mature adult")

tx_id_assignments_in = open("../../results/correctionLogs/ML_stringent.gene.txt")
ml_slope,ml_pval = plotReadCountsVsPolyA(axes[1,2],tx_id_assignments_in,read_id_to_length,"male")

tx_id_assignments_in = open("../../results/correctionLogs/all_stringent.gene.txt")
_,all_pval = plotReadCountsVsPolyA(axes[1,3],tx_id_assignments_in,read_id_to_length,"all")

_,corrected_pvals,_,_ = multipletests([l1_pval,l2_pval,l3_pval,l4_pval,ya_pval,ga_pval,ml_pval,all_pval],method="fdr_bh")

plotPval(axes[0,0],corrected_pvals[0])
plotPval(axes[0,1],corrected_pvals[1])
plotPval(axes[0,2],corrected_pvals[2])
plotPval(axes[0,3],corrected_pvals[3])
plotPval(axes[1,0],corrected_pvals[4])
plotPval(axes[1,1],corrected_pvals[5])
plotPval(axes[1,2],corrected_pvals[6])
plotPval(axes[1,3],corrected_pvals[7])

plt.tight_layout()
# plt.savefig("../../figures/supplementals/sfigure5/sfigure5B_stringent.png",dpi=450)
plt.savefig("../../figures/supplementals/sfigure5/sfigure5B.png",dpi=450)
plt.clf()

fig, axes = plt.subplots(1,2,figsize=(5.0,2),gridspec_kw = {'width_ratios':[1, 2]})
axes[0].set_ylabel("Median Poly(A) tail length")
axes[0].set_xlabel("Log$_{10}$ Read Counts")
axes[1].set_ylabel("Regression line slope")
axes[0].spines['right'].set_visible(False)
axes[0].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
axes[1].spines['top'].set_visible(False)
tx_id_assignments_in = open("../../results/correctionLogs/all_stringent.gene.txt")
plotReadCountsVsPolyALeftPanel(axes[0],tx_id_assignments_in,read_id_to_length,"")
plotSlopes(axes[1],l1_slope,l2_slope,l3_slope,l4_slope,ya_slope,ga_slope,ml_slope,regression_color)
plt.tight_layout()
# plt.savefig("../../figures/figure4/figure4C_left_and_center_stringent.png",dpi=450)
plt.savefig("../../figures/figure4/figure4C_left_and_center.png",dpi=450)
