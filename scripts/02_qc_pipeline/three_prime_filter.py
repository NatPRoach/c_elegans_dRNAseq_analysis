#!/usr/bin/env python2

import sys
import pysam
import sets
import regex as re
import matplotlib.pyplot as plt
from Bio import Seq

def splitORF(input_seq):
    startP = re.compile('ATG')
    nuc = input_seq.replace('\n','')
    longest = (0,)
    for m in startP.finditer(nuc, overlapped=True):
        mod = len(Seq.Seq(nuc)[m.start():]) % 3
        if mod != 0:
            pro = Seq.Seq(nuc)[m.start():-mod].translate(to_stop=True)
        else:
            pro = Seq.Seq(nuc)[m.start():].translate(to_stop=True)
        if len(pro) > longest[0]:
            flag = True
            stop = nuc[m.start()+len(pro)*3: m.start()+len(pro)*3 +3 ]
            if len(stop) != 3:
                flag = False
            if stop != "TAG" and stop != "TAA" and stop != "TGA":
                flag = False
            if flag:
                longest = (len(pro), 
                           m.start(), 
                           m.start()+len(pro)*3+3,
                           str(pro),
                           len(nuc[m.start()+len(pro)*3+3:]))# everything after stop codon
    return longest

infile  = open(sys.argv[1],'r')
stop_outfile = open(sys.argv[2],'w')
no_stop_outfile = open(sys.argv[3],'w')
# out_file.write("gene_id\tfull_sequence\tCDS\tCDS_start\tCDS_end\t#ofexons\n")
gene_to_isoforms = {}
gene_to_stop_codons = {}
counter = 0
### 4 - create a structure that groups read info by their associated gene
gene_to_reads = {}
for line in infile:
    fields = line.strip().split()
    longest = splitORF(fields[-1])
    if longest == (0,):
        no_stop_outfile.write(line)
    else:
        stop_outfile.write(line)