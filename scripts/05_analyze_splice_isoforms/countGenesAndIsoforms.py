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
def processNote(note):
    d = {}
    fields = note.split(" %3B")
    # print fields
    s = sets.Set()
    for field in fields:
        if field == '':
            continue
        key,val = field.strip().split()
        # if val[0:2] == "yk":
#             val = val[:-2]
            # print val
        if key in d:
            if val not in s:
                d[key].append(val)
        else:
            d[key] = [val]
        s.add(val)
    return d

cDNA_gff_in = open("/Users/nproach/Documents/NPR_Notebook/00_Data/references/ce11/c_elegans.PRJNA13758.WS265.Confirmed_cDNA.gff3")
intron_gff_in = open("/Users/nproach/Documents/NPR_Notebook/00_Data/references/ce11/c_elegans.PRJNA13758.WS265.introns.gff3")
exon_gff_in = open("/Users/nproach/Documents/NPR_Notebook/00_Data/references/ce11/c_elegans.PRJNA13758.WS265.exon.gff3")
mRNA_in = open("/Users/nproach/Documents/NPR_Notebook/00_Data/references/ce11/c_elegans.PRJNA13758.WS265.mRNA.gff3")

##For now assume an isoform has full length suport if any of its introns are supported by a cDNA (bad assumption we can readress later)

supported_intron_count = {}
all_intron_count = {}
supported_isoforms = sets.Set()
supported_genes = sets.Set()
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
for line in cDNA_gff_in:
    fields = line.strip().split('\t')
    attributes = splitAttributes(fields[8])
    tx = attributes["Parent"]
    tx_id = tx.split(':')[1]
    if tx_id in supported_intron_count:
        supported_intron_count[tx_id] += 1
    else:
        supported_intron_count[tx_id] = 1
    # uniq_isoforms.add(tx_id)
    # if tx_id in tx_id_to_gene_id:
    #     gene_id = tx_id_to_gene_id[tx_id]
    #     uniq_genes.add(gene_id)
    # else:
    #     print "something went wrong"
    #     print tx_id

# for line in intron_gff_in:
#     fields = line.strip().split('\t')
#     attributes = splitAttributes(fields[8])
#     tx = attributes["Parent"]
#     tx_id = tx.split(':')[1]
#     all_isoforms.add(tx_id)
#     if tx_id in all_intron_count:
#         all_intron_count[tx_id] += 1
#     else:
#         all_intron_count[tx_id] = 1
#     if tx_id in tx_id_to_gene_id:
#         gene_id = tx_id_to_gene_id[tx_id]
#         all_genes.add(gene_id)
#     else:
#         print "Something went wrong2"
#         print tx_id
all_intron_count = {}
intron_gff_in = open("/Users/nproach/Documents/NPR_Notebook/00_Data/references/ce11/c_elegans.PRJNA13758.WS265.introns.gff3")
for line in intron_gff_in:
    fields = line.strip().split('\t')
    attributes = splitAttributes(fields[8])
    tx = attributes["Parent"]
    tx_id = tx.split(':')[1]
    all_isoforms.add(tx_id)
    if "Note" in attributes:
        note = attributes["Note"]
        support = processNote(note)
        if tx_id not in all_intron_count:
            all_intron_count[tx_id] = [0,{}]
        if "Confirmed_EST" in support:
            for est in support["Confirmed_EST"]:
                if est in all_intron_count[tx_id][1]:
                    all_intron_count[tx_id][1][est] += 1
                else:
                    all_intron_count[tx_id][1][est] = 1
        if "Confirmed_cDNA" in support:
            for cdna in support["Confirmed_cDNA"]:
                if cdna in all_intron_count[tx_id][1]:
                    all_intron_count[tx_id][1][cdna] += 1
                else:
                    all_intron_count[tx_id][1][cdna] = 1
        # print support
    if tx_id in all_intron_count:
        all_intron_count[tx_id][0] += 1
    else:
        all_intron_count[tx_id] = [1,{}]
for tx_id in all_intron_count:
    total = all_intron_count[tx_id][0]
    for est_or_cdna in all_intron_count[tx_id][1]:
        assert all_intron_count[tx_id][1][est_or_cdna] <= total
        if all_intron_count[tx_id][1][est_or_cdna] == total:
            supported_isoforms.add(tx_id)
            if tx_id in tx_id_to_gene_id:
                gene_id = tx_id_to_gene_id[tx_id]
                supported_genes.add(gene_id)
            else:
                print "Something went wrong3"
                print tx_id
            # print "here"
            break

    if tx_id in tx_id_to_gene_id:
        gene_id = tx_id_to_gene_id[tx_id]
        all_genes.add(gene_id)
    else:
        print "Something went wrong2"
        print tx_id
# for tx_id in supported_intron_count:
#     if supported_intron_count[tx_id] == all_intron_count[tx_id]:
#         supported_isoforms.add(tx_id)
#         if tx_id in tx_id_to_gene_id:
#             gene_id = tx_id_to_gene_id[tx_id]
#             supported_genes.add(gene_id)
#         else:
#             print "Something went wrong3"
#             print tx_id
    # else:
    #     print tx_id
# print len(supported_isoforms)
    # print all_intron_count
# for line in exon_gff_in:
#     fields = line.strip().split('\t')
#     attributes = splitAttributes(fields[8])
#     tx = attributes["Parent"]
#     tx_id = tx.split(':')[1]
#     all_isoforms.add(tx_id)
#     if tx_id in tx_id_to_gene_id:
#         gene_id = tx_id_to_gene_id[tx_id]
#         all_genes.add(gene_id)
#     else:
#         print "Something went wrong2"
#         print tx_id



our_isoforms = sets.Set()
our_genes = sets.Set()
tx_id_to_introns1 = {}
tx_id_to_introns2 = {}
introns_to_tx_id = {}
intron_gff_in = open("/Users/nproach/Documents/NPR_Notebook/00_Data/references/ce11/c_elegans.PRJNA13758.WS265.introns.gff3")
for line in intron_gff_in:
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
    if fields[2] == "intron":
        attr = splitAttributes(fields[8])# .split(':')[1]
        tx_id = attr["Parent"].split(':')[1]
        # print index
        # print tx_id
        # assert tx_id in tx_id_to_gene
        if tx_id in tx_id_to_introns1:
            tx_id_to_introns1[tx_id].append((chrom,strand,start,end))
        else:
            tx_id_to_introns1[tx_id] = [(chrom,strand,start,end)]
for tx_id in tx_id_to_introns1:
    introns1 = tx_id_to_introns1[tx_id]
    introns2 = []
    for chrom,strand,start,end in introns1:
        if strand == '+':
            donor = start
            acceptor = end
        elif strand == '-':
            donor = end
            acceptor = start
        introns2.append((donor,acceptor))
    if strand == '+':
        introns2 = tuple(introns2)
    elif strand == '-':
        introns2.reverse()
        introns2 = tuple(introns2)
    if introns2 in introns_to_tx_id:
        introns_to_tx_id[introns2].append(tx_id)
    else:
        introns_to_tx_id[introns2] = [tx_id]
for introns in introns_to_tx_id:
    for tx_id in introns_to_tx_id[introns]:
        tx_id_to_introns2[tx_id] = introns

infile = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_isoforms/all_isoforms.tsv")
outfile1 = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/countGenesAndIsoforms/full_length_overlap.matrix",'w')

outfile1.write("class\ttx_id\tgene_id\tsupport\n")
uniq_introns = sets.Set()
for line in infile:
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
            assert introns not in uniq_introns
            uniq_introns.add(introns)
            if introns in introns_to_tx_id:
                tx_id = ",".join(introns_to_tx_id[introns])
                # if tx_id not in our_isoforms:
                #     print "Our isoforms\t%s\t%s\tfull length" %(tx_id,gene_id)
                outfile1.write("Our isoforms\t%s\t%s\tfull length\n" %(tx_id,gene_id))
                our_isoforms.add(tx_id)
                if gene_id not in our_genes:
                    our_genes.add(gene_id)
                    outfile1.write("Our genes\t-\t%s\tfull length\n" %(gene_id))
                
            # else:
                # print line.strip() # Report any novel isoforms by printing them to stdout
            #print introns, "reads"
            # total += 1


uniq_supported_introns = sets.Set()
uniq_unsupported_introns = sets.Set()
supported_genes = sets.Set()
# unsupported_genes = sets.Set()
for tx_id in supported_isoforms:
    introns = tx_id_to_introns2[tx_id]
    uniq_supported_introns.add(introns)
for tx_id in all_isoforms.difference(supported_isoforms):
    introns = tx_id_to_introns2[tx_id]
    uniq_unsupported_introns.add(introns)

supported_isoforms = sets.Set()
for introns in uniq_supported_introns:
    tx_id = introns_to_tx_id[introns][0]
    gene_id = tx_id_to_gene_id[tx_id]
    tx_id = ','.join(introns_to_tx_id[introns])
    supported_genes.add(gene_id)
    outfile1.write("WB isoforms\t%s\t%s\tfull length\n" %(tx_id,gene_id))
    supported_isoforms.add(tx_id)
    
for introns in uniq_unsupported_introns:
    tx_id = introns_to_tx_id[introns][0]
    gene_id = tx_id_to_gene_id[tx_id]
    tx_id = ','.join(introns_to_tx_id[introns])
    # unsupported_genes.add(gene_id)
    outfile1.write("WB isoforms\t%s\t%s\tinferred\n" %(tx_id,gene_id))
    
for gene_id in supported_genes:
    outfile1.write("WB genes\t-\t%s\tfull length\n" %(gene_id))
unsupported_genes = all_genes.difference(supported_genes)
for gene_id in unsupported_genes:
    outfile1.write("WB genes\t-\t%s\tinferred\n" %(gene_id))
# print len(supported_genes.intersection(unsupported_genes))
assert len(supported_genes.intersection(unsupported_genes)) == 0


outfile2 = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/countGenesAndIsoforms/full_length_gene_overlap.matrix",'w')
gene_union = supported_genes.union(our_genes)
for gene in gene_union:
    outfile2.write("%i\t%i\n" %(gene in our_genes, gene in supported_genes))

outfile3 = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/countGenesAndIsoforms/full_length_isoform_overlap.matrix",'w')
isoform_union = supported_isoforms.union(our_isoforms)
for isoform in isoform_union:
    outfile3.write("%i\t%i\n" %(isoform in our_isoforms, isoform in supported_isoforms))



# for tx_id in supported_isoforms:
#     gene_id = tx_id_to_gene_id[tx_id]
#     print "WormBase isoforms\t%s\t%s\tfull length" %(tx_id,gene_id)
#
# for tx_id in all_isoforms.difference(supported_isoforms):
#     gene_id = tx_id_to_gene_id[tx_id]
#     print "WormBase isoforms\t%s\t%s\tinferred" %(tx_id,gene_id)

# for gene_id in supported_genes:
#     print "WormBase genes\t-\t%s\tfull length" %(gene_id)
#
# for gene_id in all_genes.difference(supported_genes):
#     print "WormBase genes\t-\t%s\tinferred" %(gene_id)

#print previous_existing, unambiguous, total

# print len(tx_id_to_gene_id)
# print len(supported_isoforms)
# print len(all_isoforms)
# print len(supported_genes)
# print len(all_genes)
