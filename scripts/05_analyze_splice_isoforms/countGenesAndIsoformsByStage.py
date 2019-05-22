#!/usr/bin/env python2

import sets
import sys

def splitAttributes(attr):
    d = {}
    fields = attr.split(';')
    for field in fields:
        key,val = field.split('=')
        d[key] = val
    return d
def splitFastaFields(fields):
    d = {}
    for field in fields:
        key,val = field.split('=')
        d[key] = val
    return d

mRNA_in = open("../../references/WS265/c_elegans.PRJNA13758.WS265.mRNA.gff3")
outfile1 = open("../../results/scratch/countGenesAndIsoforms/staged_gene_and_isoform_count.txt",'w')
outfile2 = open("../../results/scratch/countGenesAndIsoforms/staged_novel_gene_and_isoform_count.txt",'w')

threshold = 1.0
supported_intron_count = {}
all_intron_count = {}
all_isoforms = sets.Set()
all_genes = sets.Set()

tx_id_to_gene_id = {}

for line in mRNA_in:
    fields = line.strip().split('\t')
    attributes = splitAttributes(fields[8])
    g = attributes["Parent"]
    gene_id = g.split(':')[1]
    tx = attributes["ID"]
    tx_id = tx.split(':')[1]
    tx_id_to_gene_id[tx_id] = gene_id
    all_genes.add(gene_id)



our_isoforms = sets.Set()
our_genes = sets.Set()
tx_id_to_introns1 = {}
tx_id_to_introns2 = {}
introns_to_tx_id = {}

### 1 - get annotation info, the intron chain for each transcript isoform
ingff = open("../../references/WS265/c_elegans.PRJNA13758.WS265.exon.gff3")
tx_id_to_exons = {}
for line in ingff:
    fields = line.strip().split()
    chrom = fields[0]
    if chrom == "MtDNA":
        chrom = "chrM"
    else:
        if len(chrom) < 3:
            chrom = "chr" + chrom
        elif chrom[:3] != "chr":
            chrom = "chr" + chrom
    start = int(fields[3]) - 1
    end =  int(fields[4])
    strand = fields[6]
    if fields[2] == "exon":
        tx_id = fields[8].split(':')[1]
        if tx_id in tx_id_to_exons:
            tx_id_to_exons[tx_id].append((chrom,strand,start,end))
        else:
            tx_id_to_exons[tx_id] = [(chrom,strand,start,end)]
    

tx_id_to_introns = {}
introns_to_tx_id = {}
for tx_id in tx_id_to_exons: ##WARNING: GFF3 file must be sorted for this approach to work
    exons = tx_id_to_exons[tx_id]
    if len(exons) > 1:
        introns = []
        chrom = None
        strand = None
        for x in range(len(exons)-1):
            exon1 = exons[x]
            exon2 = exons[x+1]
            assert exon1[0] == exon2[0]
            chrom = exon1[0]
            assert exon1[1] == exon2[1]
            strand = exon1[1]
            if strand == '+':
                donor = exon1[3]
                acceptor = exon2[2]
            elif strand == '-':
                donor = exon2[2]
                acceptor = exon1[3]
            introns.append((donor,acceptor))
        if strand == '+':
            introns = tuple(introns)
        elif strand == '-':
            introns.reverse()
            introns = tuple(introns)
        #store intron structure, and all 5' truncations of intron structure and link to the corresponding tx_id
        tx_id_to_introns[tx_id] = introns
        for x in range(len(introns)):
            intron_subset = introns[x:]
            if intron_subset in introns_to_tx_id:
                introns_to_tx_id[intron_subset].append(tx_id)
            else:
                introns_to_tx_id[intron_subset] = [tx_id]

l1_isoform_count = 0
l2_isoform_count = 0
l3_isoform_count = 0
l4_isoform_count = 0
ya_isoform_count = 0
ga_isoform_count = 0
ml_isoform_count = 0
all_isoform_count = 0


l1_gene_count = 0
l2_gene_count = 0
l3_gene_count = 0
l4_gene_count = 0
ya_gene_count = 0
ga_gene_count = 0
ml_gene_count = 0
all_gene_count = 0


l1_in = open("../../results/isoforms/L1_isoforms.tsv")
l2_in = open("../../results/isoforms/L2_isoforms.tsv")
l3_in = open("../../results/isoforms/L3_isoforms.tsv")
l4_in = open("../../results/isoforms/L4_isoforms.tsv")
ya_in = open("../../results/isoforms/YA_isoforms.tsv")
ga_in = open("../../results/isoforms/GA_isoforms.tsv")
ml_in = open("../../results/isoforms/ML_isoforms.tsv")
all_in = open("../../results/isoforms/all_isoforms.tsv")

l1_novel_count = 0
l2_novel_count = 0
l3_novel_count = 0
l4_novel_count = 0
ya_novel_count = 0
ga_novel_count = 0
ml_novel_count = 0
all_novel_count = 0

l1_novel_genes = sets.Set()
l2_novel_genes = sets.Set()
l3_novel_genes = sets.Set()
l4_novel_genes = sets.Set()
ya_novel_genes = sets.Set()
ga_novel_genes = sets.Set()
ml_novel_genes = sets.Set()
all_novel_genes = sets.Set()

our_genes = sets.Set()
for line in l1_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        if fields[14] != '1': # Ignore intron retained reads
            intron_string = fields[15]
            introns = []
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            if introns in introns_to_tx_id:
                tx_id = ",".join(introns_to_tx_id[introns])
                l1_isoform_count += 1
                if gene_id not in our_genes:
                    our_genes.add(gene_id)
                    l1_gene_count += 1
            else:
                l1_novel_count +=1
                l1_novel_genes.add(gene_id)
    else: #single exon gene:
        l1_isoform_count += 1
        if gene_id not in our_genes and gene_id in all_genes:
            our_genes.add(gene_id)
            l1_gene_count += 1

our_genes = sets.Set()
for line in l2_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        if fields[14] != '1': # Ignore intron retained reads
            intron_string = fields[15]
            introns = []
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            if introns in introns_to_tx_id:
                tx_id = ",".join(introns_to_tx_id[introns])
                l2_isoform_count += 1
                if gene_id not in our_genes:
                    our_genes.add(gene_id)
                    l2_gene_count += 1
            else:
                l2_novel_count +=1
                l2_novel_genes.add(gene_id)
    else: #single exon gene:
        l2_isoform_count += 1
        if gene_id not in our_genes and gene_id in all_genes:
            our_genes.add(gene_id)
            l2_gene_count += 1

our_genes = sets.Set()
for line in l3_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        if fields[14] != '1': # Ignore intron retained reads
            intron_string = fields[15]
            introns = []
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            if introns in introns_to_tx_id:
                tx_id = ",".join(introns_to_tx_id[introns])
                l3_isoform_count += 1
                if gene_id not in our_genes:
                    our_genes.add(gene_id)
                    l3_gene_count += 1
            else:
                l3_novel_count +=1
                l3_novel_genes.add(gene_id)
    else: #single exon gene:
        l3_isoform_count += 1
        if gene_id not in our_genes and gene_id in all_genes:
            our_genes.add(gene_id)
            l3_gene_count += 1

                
our_genes = sets.Set()
for line in l4_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        if fields[14] != '1': # Ignore intron retained reads
            intron_string = fields[15]
            introns = []
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            if introns in introns_to_tx_id:
                tx_id = ",".join(introns_to_tx_id[introns])
                l4_isoform_count += 1
                if gene_id not in our_genes:
                    our_genes.add(gene_id)
                    l4_gene_count += 1
            else:
                l4_novel_count +=1
                l4_novel_genes.add(gene_id)
    else: #single exon gene:
        l4_isoform_count += 1
        if gene_id not in our_genes and gene_id in all_genes:
            our_genes.add(gene_id)
            l4_gene_count += 1
    
our_genes = sets.Set()
for line in ya_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        if fields[14] != '1': # Ignore intron retained reads
            intron_string = fields[15]
            introns = []
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            if introns in introns_to_tx_id:
                tx_id = ",".join(introns_to_tx_id[introns])
                ya_isoform_count += 1
                if gene_id not in our_genes:
                    our_genes.add(gene_id)
                    ya_gene_count += 1
            else:
                ya_novel_count +=1
                ya_novel_genes.add(gene_id)
    else: #single exon gene:
        ya_isoform_count += 1
        if gene_id not in our_genes and gene_id in all_genes:
            our_genes.add(gene_id)
            ya_gene_count += 1
                    
our_genes = sets.Set()
for line in ga_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        if fields[14] != '1': # Ignore intron retained reads
            intron_string = fields[15]
            introns = []
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            if introns in introns_to_tx_id:
                tx_id = ",".join(introns_to_tx_id[introns])
                ga_isoform_count += 1
                if gene_id not in our_genes:
                    our_genes.add(gene_id)
                    ga_gene_count += 1
            else:
                ga_novel_count +=1
                ga_novel_genes.add(gene_id)
    else: #single exon gene:
        ga_isoform_count += 1
        if gene_id not in our_genes and gene_id in all_genes:
            our_genes.add(gene_id)
            ga_gene_count += 1
        
our_genes = sets.Set()
for line in ml_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        if fields[14] != '1': # Ignore intron retained reads
            intron_string = fields[15]
            introns = []
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            if introns in introns_to_tx_id:
                tx_id = ",".join(introns_to_tx_id[introns])
                ml_isoform_count += 1
                if gene_id not in our_genes:
                    our_genes.add(gene_id)
                    ml_gene_count += 1
            else:
                ml_novel_count +=1
                ml_novel_genes.add(gene_id)
    else: #single exon gene:
        ml_isoform_count += 1
        if gene_id not in our_genes and gene_id in all_genes:
            our_genes.add(gene_id)
            ml_gene_count += 1

our_genes = sets.Set()
for line in all_in:
    fields = line.strip().split()
    gene_id = fields[0]
    if len(fields) == 16:
        if fields[14] != '1': # Ignore intron retained reads
            intron_string = fields[15]
            introns = []
            for pair in intron_string.split(';'):
                split_pair = pair.split(',')
                if len(split_pair) == 2:
                    intron = (int(split_pair[0]),int(split_pair[1]))
                    introns.append(intron)
            introns = tuple(introns)
            if introns in introns_to_tx_id:
                tx_id = ",".join(introns_to_tx_id[introns])
                all_isoform_count += 1
                if gene_id not in our_genes:
                    our_genes.add(gene_id)
                    all_gene_count += 1
            else:
                all_novel_count +=1
                all_novel_genes.add(gene_id)
    else: #single exon gene:
        all_isoform_count += 1
        if gene_id not in our_genes and gene_id in all_genes:
            our_genes.add(gene_id)
            all_gene_count += 1


outfile1.write("stage\tdataset\tcounts\tx_order\n")
outfile1.write("L1\tOur isoforms\t%d\t0\n" %(l1_isoform_count))
outfile1.write("L2\tOur isoforms\t%d\t1\n" %(l2_isoform_count))
outfile1.write("L3\tOur isoforms\t%d\t2\n" %(l3_isoform_count))
outfile1.write("L4\tOur isoforms\t%d\t3\n" %(l4_isoform_count))
outfile1.write("young adult\tOur isoforms\t%d\t4\n" %(ya_isoform_count))
outfile1.write("mature adult\tOur isoforms\t%d\t5\n" %(ga_isoform_count))
outfile1.write("male\tOur isoforms\t%d\t6\n" %(ml_isoform_count))
outfile1.write("all\tOur isoforms\t%d\t7\n" %(all_isoform_count))

outfile1.write("L1\tOur genes\t%d\t0\n" %(l1_gene_count))
outfile1.write("L2\tOur genes\t%d\t1\n" %(l2_gene_count))
outfile1.write("L3\tOur genes\t%d\t2\n" %(l3_gene_count))
outfile1.write("L4\tOur genes\t%d\t3\n" %(l4_gene_count))
outfile1.write("young adult\tOur genes\t%d\t4\n" %(ya_gene_count))
outfile1.write("mature adult\tOur genes\t%d\t5\n" %(ga_gene_count))
outfile1.write("male\tOur genes\t%d\t6\n" %(ml_gene_count))
outfile1.write("all\tOur genes\t%d\t7\n" %(all_gene_count))

outfile2.write("stage\tdataset\tcounts\tx_order\n")
outfile2.write("L1\tnovel isoforms\t%s\t0\n" %(l1_novel_count))
outfile2.write("L2\tnovel isoforms\t%s\t1\n" %(l2_novel_count))
outfile2.write("L3\tnovel isoforms\t%s\t2\n" %(l3_novel_count))
outfile2.write("L4\tnovel isoforms\t%s\t3\n" %(l4_novel_count))
outfile2.write("young adult\tnovel isoforms\t%s\t4\n" %(ya_novel_count))
outfile2.write("mature adult\tnovel isoforms\t%s\t5\n" %(ga_novel_count))
outfile2.write("male\tnovel isoforms\t%s\t6\n" %(ml_novel_count))
outfile2.write("all\tnovel isoforms\t%s\t7\n" %(all_novel_count))

outfile2.write("L1\tgenes with novel isoforms\t%s\t0\n" %(len(l1_novel_genes)))
outfile2.write("L2\tgenes with novel isoforms\t%s\t1\n" %(len(l2_novel_genes)))
outfile2.write("L3\tgenes with novel isoforms\t%s\t2\n" %(len(l3_novel_genes)))
outfile2.write("L4\tgenes with novel isoforms\t%s\t3\n" %(len(l4_novel_genes)))
outfile2.write("young adult\tgenes with novel isoforms\t%s\t4\n" %(len(ya_novel_genes)))
outfile2.write("mature adult\tgenes with novel isoforms\t%s\t5\n" %(len(ga_novel_genes)))
outfile2.write("male\tgenes with novel isoforms\t%s\t6\n" %(len(ml_novel_genes)))
outfile2.write("all\tgenes with novel isoforms\t%s\t7\n" %(len(all_novel_genes)))