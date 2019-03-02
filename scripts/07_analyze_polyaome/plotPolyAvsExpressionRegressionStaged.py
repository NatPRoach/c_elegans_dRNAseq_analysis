#!/usr/bin/env python2
import sys
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
#from sklearn.pipeline import Pipeline

def load_read_id_to_length():
    read_id_to_length = {}### For each read in the polyA file assign its polyA tail length in a dictionary
    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L1/ce11_gen_L1_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L1/ce11_gen_L1_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L2/ce11_gen_L2_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L2/ce11_gen_L2_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])            

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L3/ce11_gen_L3_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L3/ce11_gen_L3_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L4/ce11_gen_L4_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/L4/ce11_gen_L4_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])
            
    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/YA/ce11_gen_young_adult_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/YA/ce11_gen_young_adult_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/GA/ce11_gen_adult_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/GA/ce11_gen_adult_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])
            
    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/ML/ce11_gen_male_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("/Users/nproach/Documents/NPR_Notebook/00_Data/c_elegans/analysis/cross_stage/polya/ML/ce11_gen_male_bio1_tech2_Minimap2_k14.polya",'r')
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
        if len(polyas) != 0:
            median_polya.append(np.median(polyas))
            log_expression.append(np.log2(float(len(read_ids))))

    matrix = np.matrix([log_expression,median_polya])
    matrix = matrix.transpose()
    # poly = PolynomialFeatures(degree=2)
    # model = Pipeline([('poly', PolynomialFeatures(degree=2)),('linear', LinearRegression(fit_intercept=False))])
    model = LinearRegression()
    model = model.fit(matrix[:,0],matrix[:,1])
    coefs = [model.intercept_[0],model.coef_[0][0]]
    print coefs
    #coefs = model.named_steps['linear'].coef_
    x = np.linspace(0,12,1000)
    y = coefs[0] + coefs[1]*x #+ coefs[:,2]*(x**2)
    y_hat = coefs[0] + coefs[1]*np.array(log_expression)
    r2 = r2_score(np.array(median_polya),y_hat)
    print r2
    ax.scatter(log_expression,median_polya,alpha = 0.05,color=point_color)
    ax.plot(x,y,color=regression_color)
    ax.text(6,150,"$R^2$=%.4f"%(r2))
    # ax.ylabel("Median PolyA tail length")
    # ax.xlabel("Log2 Read Counts")
    ax.set_ylim(0,200)
    ax.set_title(title)
    
    
    # plt.scatter(log_expression,median_polya,alpha = 0.05,color="#2579B2")
#     plt.plot(x,y,color="#389E34")
#     plt.text(8,150,"$R^2$=%f"%(r2))
#     plt.ylabel("Median PolyA tail length")
#     plt.xlabel("Log2 Read Counts")
#     plt.ylim(0,200)
    #plt.show()

font = {"family":"sans-serif",
        "sans-serif":["Helvetica"],
        "size":8}
matplotlib.rc('font',**font)

if len(sys.argv) > 1:
    point_color = sys.argv[1]
    regression_color= sys.argv[2]
else:
    point_color = "#2579B2"
    regression_color = "#FD7F23"

read_id_to_length = load_read_id_to_length()
#plt.figure(figsize=(20,8))
fig, axes = plt.subplots(2,4,sharex="col",sharey="row",figsize=(8.333,5.375))
axes[0,0].set_ylabel("Median PolyA tail length")
axes[1,0].set_ylabel("Median PolyA tail length")
axes[1,0].set_xlabel("Log2 Read Counts")
axes[1,1].set_xlabel("Log2 Read Counts")
axes[1,2].set_xlabel("Log2 Read Counts")
axes[1,3].set_xlabel("Log2 Read Counts")

tx_id_assignments_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/correctionLogs/L1.gene.txt")
plotReadCountsVsPolyA(axes[0,0],tx_id_assignments_in,read_id_to_length,"L1")

tx_id_assignments_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/correctionLogs/L2.gene.txt")
plotReadCountsVsPolyA(axes[0,1],tx_id_assignments_in,read_id_to_length,"L2")

tx_id_assignments_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/correctionLogs/L3.gene.txt")
plotReadCountsVsPolyA(axes[0,2],tx_id_assignments_in,read_id_to_length,"L3")

tx_id_assignments_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/correctionLogs/L4.gene.txt")
plotReadCountsVsPolyA(axes[0,3],tx_id_assignments_in,read_id_to_length,"L4")

tx_id_assignments_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/correctionLogs/YA.gene.txt")
plotReadCountsVsPolyA(axes[1,0],tx_id_assignments_in,read_id_to_length,"young adult")

tx_id_assignments_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/correctionLogs/GA.gene.txt")
plotReadCountsVsPolyA(axes[1,1],tx_id_assignments_in,read_id_to_length,"mature adult")

tx_id_assignments_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/correctionLogs/ML.gene.txt")
plotReadCountsVsPolyA(axes[1,2],tx_id_assignments_in,read_id_to_length,"male")
# plt.ylabel("Median PolyA tail length")
# plt.xlabel("Log2 Read Counts")
#plt.ylim(0,200)
tx_id_assignments_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/correctionLogs/all.gene.txt")
plotReadCountsVsPolyA(axes[1,3],tx_id_assignments_in,read_id_to_length,"all")
plt.tight_layout()
#plt.savefig("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/polya_correlations/expression/plots/polyAvsExpression.pdf")

plt.savefig("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/polya_correlations/expression/plots/polyAvsExpression.png",dpi=450)
#plt.show()