#!/usr/bin/env python

infile = open("../../results/coding_predictions/all_stringent_isoforms_predictions.tsv")
outfile = open("../../results/coding_predictions/Supplemental_Table_S5.csv",'w')

outfile.write("Spliceform ID,\n")
# max_reads = 0
# for i,line in enumerate(infile):
#     if i == 0:
#         continue
#     fields = line.strip().split()
#     if float(fields[-1]) < 0.5:
#         gene,read_ids = fields[0].split('_')
#         max_reads = max((len(read_ids.split(',')), max_reads))
#
# infile = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/isoforms/all_stringent_isoforms.tsv")


count = 0
gene_spliceform_counter = {}
for i,line in enumerate(infile):
    if i == 0:
        continue
    fields = line.strip().split()
    if float(fields[-1]) < 0.5:
        gene,read_ids = fields[0].split('_')
        gene = "WBGene" + gene[6:]
        if gene not in gene_spliceform_counter:
            gene_spliceform_counter[gene] = 0
        else:
            gene_spliceform_counter[gene] += 1
        count += 1
        outfile.write("%s,\n" %(gene + "-Spliceform" + str(gene_spliceform_counter[gene])))
    # else:
    #     gene,read_ids = fields[0].split('_')
print count