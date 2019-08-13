#!/usr/bin/env python2

import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns
import sets


if len(sys.argv) > 1:
    color = sys.argv[1]
    color2 =sys.argv[2]
else:
    color = "#6A4098"
    color2= "#2579B2"
font = {"family":"sans-serif",
        "sans-serif":["Helvetica"],
        "size":10}
matplotlib.rc('font',**font)

def splitAttributes(attr):
    fields = attr.split(';')
    attr_dict = {}
    for field in fields:
        if len(field) == 0:
            continue
        key,val = field.strip().split(' ')
        val = val.strip('\"')
        attr_dict[key] = val
    return attr_dict

infile = open("../../references/WBcel235/Caenorhabditis_elegans.WBcel235.93.gtf")

gene_to_transcripts = {}
for line in infile:
    if line[0:2] == '#!':
        continue
    fields = line.strip().split('\t')
    chrom = fields[0]
    if chrom == "MtDNA":
        chrom = "chrM"
    else:
        chrom = "chr" + chrom
    start = int(fields[3]) - 1
    end =  int(fields[4])
    strand = fields[6]
    if fields[2] == "gene":
        attr = splitAttributes(fields[8])
        if "gene_id" in attr:
            gene_id = attr["gene_id"]
        else:
            print "something went wrong"
    if fields[2] == "transcript":
        attr = splitAttributes(fields[8])
        if "gene_id" in attr:
            gene_id = attr["gene_id"]
        else:
            print "something went wrong"
        if "transcript_id" in attr:
            tx_id = attr["transcript_id"]
        else:
            print "something went wrong"
        if "transcript_biotype" in attr:
            biotype = attr["transcript_biotype"]
        else:
            print "something went wrong"
        if biotype == "piRNA":
            continue
        if biotype == "miRNA":
            continue
        if biotype == "tRNA":
            continue
        if biotype == "snoRNA":
            continue
        if biotype == "rRNA":
            continue
        if biotype == "pseudogene":
            continue
        if biotype == "ncRNA":
            continue
        if gene_id in gene_to_transcripts:
            gene_to_transcripts[gene_id].add(tx_id)
        else:
            gene_to_transcripts[gene_id] = sets.Set([tx_id])
        

isoform_counts = []
print len(gene_to_transcripts)
for gene_id in gene_to_transcripts:
    isoform_counts.append(len(gene_to_transcripts[gene_id]))


# l1_in = open("../../results/isoforms/L1_isoforms.tsv")
# l2_in = open("../../results/isoforms/L2_isoforms.tsv")
# l3_in = open("../../results/isoforms/L3_isoforms.tsv")
# l4_in = open("../../results/isoforms/L4_isoforms.tsv")
# ya_in = open("../../results/isoforms/YA_isoforms.tsv")
# ga_in = open("../../results/isoforms/GA_isoforms.tsv")
# ml_in = open("../../results/isoforms/ML_isoforms.tsv")
# all_in = open("../../results/isoforms/all_isoforms.tsv")
#
# l1_genes = {}
# l2_genes = {}
# l3_genes = {}
# l4_genes = {}
# ya_genes = {}
# ga_genes = {}
# ml_genes = {}
# all_genes = {}
#
# for line in l1_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in l1_genes:
#         l1_genes[gene_id] += 1
#     else:
#         l1_genes[gene_id] = 1
#
# for line in l2_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in l2_genes:
#         l2_genes[gene_id] += 1
#     else:
#         l2_genes[gene_id] = 1
#
# for line in l3_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in l3_genes:
#         l3_genes[gene_id] += 1
#     else:
#         l3_genes[gene_id] = 1
#
# for line in l4_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in l4_genes:
#         l4_genes[gene_id] += 1
#     else:
#         l4_genes[gene_id] = 1
#
# for line in ya_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in ya_genes:
#         ya_genes[gene_id] += 1
#     else:
#         ya_genes[gene_id] = 1
#
# for line in ga_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in ga_genes:
#         ga_genes[gene_id] += 1
#     else:
#         ga_genes[gene_id] = 1
#
# for line in ml_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in ml_genes:
#         ml_genes[gene_id] += 1
#     else:
#         ml_genes[gene_id] = 1
#
# for line in all_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in all_genes:
#         all_genes[gene_id] += 1
#     else:
#         all_genes[gene_id] = 1
#
# l1_counts = []
# l2_counts = []
# l3_counts = []
# l4_counts = []
# ya_counts = []
# ga_counts = []
# ml_counts = []
# all_counts = []
#
# for gene in l1_genes:
#     l1_counts.append(l1_genes[gene])
#
# for gene in l2_genes:
#     l2_counts.append(l2_genes[gene])
#
# for gene in l3_genes:
#     l3_counts.append(l3_genes[gene])
#
# for gene in l4_genes:
#     l4_counts.append(l4_genes[gene])
#
# for gene in ya_genes:
#     ya_counts.append(ya_genes[gene])
#
# for gene in ga_genes:
#     ga_counts.append(ga_genes[gene])
#
# for gene in ml_genes:
#     ml_counts.append(ml_genes[gene])
#
# for gene in all_genes:
#     all_counts.append(all_genes[gene])
#
# ### Uncomment below to make additional figures breaking down isoforms per gene by stage.
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(l1_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("L1")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/01_l1isoformsPerGene.pdf")
# # plt.clf()
# # # plt.show()
# #
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(l2_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("L2")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/02_l2isoformsPerGene.pdf")
# # plt.clf()
# #
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(l3_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("L3")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/03_l3isoformsPerGene.pdf")
# # plt.clf()
# #
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(l4_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("L4")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/04_l4isoformsPerGene.pdf")
# # plt.clf()
# #
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(ya_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("young adult")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/05_yaisoformsPerGene.pdf")
# # plt.clf()
# #
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(ga_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("mature adult")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/06_gaisoformsPerGene.pdf")
# # plt.clf()
# #
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(ml_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("male")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/07_mlisoformsPerGene.pdf")
# # plt.clf()
# #
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(all_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("all")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/08_allisoformsPerGeneNoTitle.pdf")
# # plt.clf()
#
# plt.figure(num=None,figsize=(3.5,2))
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.hist(all_counts,bins= [x - 0.5 for x in range(min(all_counts),8)],density=True,rwidth=0.6,color=color,edgecolor='black')
# #plt.title("all")
# plt.xlabel("splice isoforms per gene")
# plt.ylabel("density")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/figure2/figure2D.pdf")
# plt.clf()

# ### Sensitive
#
# l1_in = open("../../results/isoforms/L1_sensitive_isoforms.tsv")
# l2_in = open("../../results/isoforms/L2_sensitive_isoforms.tsv")
# l3_in = open("../../results/isoforms/L3_sensitive_isoforms.tsv")
# l4_in = open("../../results/isoforms/L4_sensitive_isoforms.tsv")
# ya_in = open("../../results/isoforms/YA_sensitive_isoforms.tsv")
# ga_in = open("../../results/isoforms/GA_sensitive_isoforms.tsv")
# ml_in = open("../../results/isoforms/ML_sensitive_isoforms.tsv")
# all_in = open("../../results/isoforms/all_sensitive_isoforms.tsv")
#
# l1_genes = {}
# l2_genes = {}
# l3_genes = {}
# l4_genes = {}
# ya_genes = {}
# ga_genes = {}
# ml_genes = {}
# all_genes = {}
#
# for line in l1_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in l1_genes:
#         l1_genes[gene_id] += 1
#     else:
#         l1_genes[gene_id] = 1
#
# for line in l2_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in l2_genes:
#         l2_genes[gene_id] += 1
#     else:
#         l2_genes[gene_id] = 1
#
# for line in l3_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in l3_genes:
#         l3_genes[gene_id] += 1
#     else:
#         l3_genes[gene_id] = 1
#
# for line in l4_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in l4_genes:
#         l4_genes[gene_id] += 1
#     else:
#         l4_genes[gene_id] = 1
#
# for line in ya_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in ya_genes:
#         ya_genes[gene_id] += 1
#     else:
#         ya_genes[gene_id] = 1
#
# for line in ga_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in ga_genes:
#         ga_genes[gene_id] += 1
#     else:
#         ga_genes[gene_id] = 1
#
# for line in ml_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in ml_genes:
#         ml_genes[gene_id] += 1
#     else:
#         ml_genes[gene_id] = 1
#
# for line in all_in:
#     fields = line.strip().split()
#     gene_id = fields[0]
#     if gene_id in all_genes:
#         all_genes[gene_id] += 1
#     else:
#         all_genes[gene_id] = 1
#
# l1_counts = []
# l2_counts = []
# l3_counts = []
# l4_counts = []
# ya_counts = []
# ga_counts = []
# ml_counts = []
# all_counts = []
#
# for gene in l1_genes:
#     l1_counts.append(l1_genes[gene])
#
# for gene in l2_genes:
#     l2_counts.append(l2_genes[gene])
#
# for gene in l3_genes:
#     l3_counts.append(l3_genes[gene])
#
# for gene in l4_genes:
#     l4_counts.append(l4_genes[gene])
#
# for gene in ya_genes:
#     ya_counts.append(ya_genes[gene])
#
# for gene in ga_genes:
#     ga_counts.append(ga_genes[gene])
#
# for gene in ml_genes:
#     ml_counts.append(ml_genes[gene])
#
# for gene in all_genes:
#     all_counts.append(all_genes[gene])
#
# ### Uncomment below to make additional figures breaking down isoforms per gene by stage.
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(l1_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("L1")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/01_l1isoformsPerGene.pdf")
# # plt.clf()
# # # plt.show()
# #
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(l2_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("L2")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/02_l2isoformsPerGene.pdf")
# # plt.clf()
# #
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(l3_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("L3")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/03_l3isoformsPerGene.pdf")
# # plt.clf()
# #
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(l4_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("L4")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/04_l4isoformsPerGene.pdf")
# # plt.clf()
# #
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(ya_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("young adult")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/05_yaisoformsPerGene.pdf")
# # plt.clf()
# #
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(ga_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("mature adult")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/06_gaisoformsPerGene.pdf")
# # plt.clf()
# #
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(ml_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("male")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/07_mlisoformsPerGene.pdf")
# # plt.clf()
# #
# # plt.figure(num=None,figsize=(3.5,2))
# # plt.hist(all_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# # plt.title("all")
# # plt.xlabel("splice isoforms per gene")
# # plt.ylabel("density")
# # plt.ylim((0,0.95))
# # plt.tight_layout()
# # plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/08_allisoformsPerGeneNoTitle.pdf")
# # plt.clf()
#
# plt.figure(num=None,figsize=(4.1667,2.5))
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.hist(all_counts,bins= [x - 0.5 for x in range(min(all_counts),8)],density=True,rwidth=0.6,color=color,edgecolor='black')
# #plt.title("all")
# plt.xlabel("splice isoforms per gene")
# plt.ylabel("proportion of genes")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/figure2/sensitive_figure2D.pdf")
# plt.clf()
#
# plt.figure(num=None,figsize=(4.1667,2.5))
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.hist([all_counts,isoform_counts],bins= [x - 0.5 for x in range(min(all_counts),8)],density=True,rwidth=0.6,color=[color,color2],edgecolor='black',label=["This study","annotation"])
# #plt.title("all")
# plt.legend(frameon=False)
# plt.xlabel("splice isoforms per gene")
# plt.ylabel("proportion of genes")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/figure2/sensitive_figure2D_modified.pdf")
# plt.clf()
#






### Stringent

l1_in = open("../../results/isoforms/L1_stringent_isoforms.tsv")
l2_in = open("../../results/isoforms/L2_stringent_isoforms.tsv")
l3_in = open("../../results/isoforms/L3_stringent_isoforms.tsv")
l4_in = open("../../results/isoforms/L4_stringent_isoforms.tsv")
ya_in = open("../../results/isoforms/YA_stringent_isoforms.tsv")
ga_in = open("../../results/isoforms/GA_stringent_isoforms.tsv")
ml_in = open("../../results/isoforms/ML_stringent_isoforms.tsv")
all_in = open("../../results/isoforms/all_stringent_isoforms.tsv")

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

### Uncomment below to make additional figures breaking down isoforms per gene by stage.
# plt.figure(num=None,figsize=(3.5,2))
# plt.hist(l1_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# plt.title("L1")
# plt.xlabel("splice isoforms per gene")
# plt.ylabel("density")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/01_l1isoformsPerGene.pdf")
# plt.clf()
# # plt.show()
#
# plt.figure(num=None,figsize=(3.5,2))
# plt.hist(l2_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# plt.title("L2")
# plt.xlabel("splice isoforms per gene")
# plt.ylabel("density")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/02_l2isoformsPerGene.pdf")
# plt.clf()
#
# plt.figure(num=None,figsize=(3.5,2))
# plt.hist(l3_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# plt.title("L3")
# plt.xlabel("splice isoforms per gene")
# plt.ylabel("density")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/03_l3isoformsPerGene.pdf")
# plt.clf()
#
# plt.figure(num=None,figsize=(3.5,2))
# plt.hist(l4_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# plt.title("L4")
# plt.xlabel("splice isoforms per gene")
# plt.ylabel("density")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/04_l4isoformsPerGene.pdf")
# plt.clf()
#
# plt.figure(num=None,figsize=(3.5,2))
# plt.hist(ya_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# plt.title("young adult")
# plt.xlabel("splice isoforms per gene")
# plt.ylabel("density")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/05_yaisoformsPerGene.pdf")
# plt.clf()
#
# plt.figure(num=None,figsize=(3.5,2))
# plt.hist(ga_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# plt.title("mature adult")
# plt.xlabel("splice isoforms per gene")
# plt.ylabel("density")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/06_gaisoformsPerGene.pdf")
# plt.clf()
#
# plt.figure(num=None,figsize=(3.5,2))
# plt.hist(ml_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# plt.title("male")
# plt.xlabel("splice isoforms per gene")
# plt.ylabel("density")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/07_mlisoformsPerGene.pdf")
# plt.clf()
#
# plt.figure(num=None,figsize=(3.5,2))
# plt.hist(all_counts,bins= [x - 0.5 for x in range(min(all_counts),11)],density=True,rwidth=0.6,color=color,edgecolor='black')
# plt.title("all")
# plt.xlabel("splice isoforms per gene")
# plt.ylabel("density")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/supplementals/isoformsPerGeneByStage/08_allisoformsPerGeneNoTitle.pdf")
# plt.clf()

# plt.figure(num=None,figsize=(4.1667,2.5))
# ax = plt.gca()
# ax.spines['right'].set_visible(False)
# ax.spines['top'].set_visible(False)
# plt.hist(all_counts,bins= [x - 0.5 for x in range(min(all_counts),8)],density=True,rwidth=0.6,color=color,edgecolor='black')
# #plt.title("all")
# plt.xlabel("splice isoforms per gene")
# plt.ylabel("proportion of genes")
# plt.ylim((0,0.95))
# plt.tight_layout()
# plt.savefig("../../figures/figure2/stringent_figure2D.pdf")
# plt.clf()
#


plt.figure(num=None,figsize=(4.1667,2.5))
ax = plt.gca()
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
plt.hist([all_counts,isoform_counts],bins= [x - 0.5 for x in range(min(all_counts),8)],density=True,rwidth=0.6,color=[color,color2],edgecolor='black',label=["This study","annotation"])
#plt.title("all")
plt.legend(frameon=False)
plt.xlabel("splice isoforms per gene")
plt.ylabel("proportion of genes")
plt.ylim((0,0.95))
plt.tight_layout()
plt.savefig("../../figures/figure2/figure2D.pdf")
plt.clf()