#!/usr/bin/env python2

import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

if len(sys.argv) > 1:
    color = sys.argv[1]
else:
    color = "#6A4098"
font = {"family":"sans-serif",
        "sans-serif":["Helvetica"],
        "size":8}
matplotlib.rc('font',**font)

l1_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_isoforms/L1_isoforms.tsv")
l2_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_isoforms/L2_isoforms.tsv")
l3_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_isoforms/L3_isoforms.tsv")
l4_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_isoforms/L4_isoforms.tsv")
ya_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_isoforms/YA_isoforms.tsv")
ga_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_isoforms/GA_isoforms.tsv")
ml_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_isoforms/ML_isoforms.tsv")
all_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_isoforms/all_isoforms.tsv")

l1_genes = {}
l2_genes = {}
l3_genes = {}
l4_genes = {}
ya_genes = {}
ga_genes = {}
ml_genes = {}
all_genes = {}

for line in l1_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if gene_id in l1_genes:
        l1_genes[gene_id] += 1
    else:
        l1_genes[gene_id] = 1

for line in l2_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if gene_id in l2_genes:
        l2_genes[gene_id] += 1
    else:
        l2_genes[gene_id] = 1

for line in l3_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if gene_id in l3_genes:
        l3_genes[gene_id] += 1
    else:
        l3_genes[gene_id] = 1

for line in l4_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if gene_id in l4_genes:
        l4_genes[gene_id] += 1
    else:
        l4_genes[gene_id] = 1
        
for line in ya_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if gene_id in ya_genes:
        ya_genes[gene_id] += 1
    else:
        ya_genes[gene_id] = 1
        
for line in ga_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if gene_id in ga_genes:
        ga_genes[gene_id] += 1
    else:
        ga_genes[gene_id] = 1
        
for line in ml_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if gene_id in ml_genes:
        ml_genes[gene_id] += 1
    else:
        ml_genes[gene_id] = 1
        
for line in all_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if gene_id in all_genes:
        all_genes[gene_id] += 1
    else:
        all_genes[gene_id] = 1

l1_counts = []
l2_counts = []
l3_counts = []
l4_counts = []
ya_counts = []
ga_counts = []
ml_counts = []
all_counts = []

for gene in l1_genes:
    l1_counts.append(l1_genes[gene])

for gene in l2_genes:
    l2_counts.append(l2_genes[gene])

for gene in l3_genes:
    l3_counts.append(l3_genes[gene])

for gene in l4_genes:
    l4_counts.append(l4_genes[gene])

for gene in ya_genes:
    ya_counts.append(ya_genes[gene])

for gene in ga_genes:
    ga_counts.append(ga_genes[gene])

for gene in ml_genes:
    ml_counts.append(ml_genes[gene])

for gene in all_genes:
    all_counts.append(all_genes[gene])

plt.figure(num=None,figsize=(3.5,2))
plt.hist(l1_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color)
plt.title("L1")
plt.xlabel("splice isoforms per gene")
plt.ylabel("density")
plt.ylim((0,0.95))
plt.tight_layout()
plt.savefig("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/countGenesAndIsoforms/plots/isoformsPerGene/01_l1isoformsPerGene.pdf")
plt.clf()
# plt.show()

plt.figure(num=None,figsize=(3.5,2))
plt.hist(l2_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color)
plt.title("L2")
plt.xlabel("splice isoforms per gene")
plt.ylabel("density")
plt.ylim((0,0.95))
plt.tight_layout()
plt.savefig("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/countGenesAndIsoforms/plots/isoformsPerGene/02_l2isoformsPerGene.pdf")
plt.clf()

plt.figure(num=None,figsize=(3.5,2))
plt.hist(l3_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color)
plt.title("L3")
plt.xlabel("splice isoforms per gene")
plt.ylabel("density")
plt.ylim((0,0.95))
plt.tight_layout()
plt.savefig("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/countGenesAndIsoforms/plots/isoformsPerGene/03_l3isoformsPerGene.pdf")
plt.clf()

plt.figure(num=None,figsize=(3.5,2))
plt.hist(l4_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color)
plt.title("L4")
plt.xlabel("splice isoforms per gene")
plt.ylabel("density")
plt.ylim((0,0.95))
plt.tight_layout()
plt.savefig("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/countGenesAndIsoforms/plots/isoformsPerGene/04_l4isoformsPerGene.pdf")
plt.clf()

plt.figure(num=None,figsize=(3.5,2))
plt.hist(ya_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color)
plt.title("young adult")
plt.xlabel("splice isoforms per gene")
plt.ylabel("density")
plt.ylim((0,0.95))
plt.tight_layout()
plt.savefig("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/countGenesAndIsoforms/plots/isoformsPerGene/05_yaisoformsPerGene.pdf")
plt.clf()

plt.figure(num=None,figsize=(3.5,2))
plt.hist(ga_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color)
plt.title("mature adult")
plt.xlabel("splice isoforms per gene")
plt.ylabel("density")
plt.ylim((0,0.95))
plt.tight_layout()
plt.savefig("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/countGenesAndIsoforms/plots/isoformsPerGene/06_gaisoformsPerGene.pdf")
plt.clf()

plt.figure(num=None,figsize=(3.5,2))
plt.hist(ml_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color)
plt.title("male")
plt.xlabel("splice isoforms per gene")
plt.ylabel("density")
plt.ylim((0,0.95))
plt.tight_layout()
plt.savefig("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/countGenesAndIsoforms/plots/isoformsPerGene/07_mlisoformsPerGene.pdf")
plt.clf()

plt.figure(num=None,figsize=(3.5,2))
plt.hist(all_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color)
plt.title("all")
plt.xlabel("splice isoforms per gene")
plt.ylabel("density")
plt.ylim((0,0.95))
plt.tight_layout()
plt.savefig("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/countGenesAndIsoforms/plots/isoformsPerGene/08_allisoformsPerGeneNoTitle.pdf")
plt.clf()

plt.figure(num=None,figsize=(3.5,2))
plt.hist(all_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color)
#plt.title("all")
plt.xlabel("splice isoforms per gene")
plt.ylabel("density")
plt.ylim((0,0.95))
plt.tight_layout()
plt.savefig("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/countGenesAndIsoforms/plots/isoformsPerGene/09_allisoformsPerGeneNoTitle.pdf")
plt.clf()