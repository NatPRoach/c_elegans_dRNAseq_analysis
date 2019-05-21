#!/usr/bin/env python2

import sys
import numpy
import matplotlib.pyplot as pyplot

def getPASBreakdown(fastaFile,outfile,altPAStable,stage,x_order):
    for line in fastaFile:
        if line[0] == '>':
            cluster = line.strip().strip(">()+-")
            #print cluster
            continue
        seq = line.strip().upper()
        canonFlag = False
        offset = 0
        for x in reversed(range(20,window-5)):
            subseq = seq[x:x+6]
            if subseq == "AATAAA":
                canonFlag = True
                break
        if canonFlag:
            outfile.write("%s\t%s\t%s\t%d\n" %(stage,cluster,"AAUAAA",x_order))
        else:
            noncanonFlag = False
            for altPAS in altPAStable:
                for x in reversed(range(20,window-5)):
                    subseq = seq[x:x+6]
                    if subseq == altPAS:
                        noncanonFlag = True
                if noncanonFlag:
                    break
            if noncanonFlag:
                outfile.write("%s\t%s\t%s\t%d\n" %(stage,cluster,"altPAS",x_order))
            else:
                outfile.write("%s\t%s\t%s\t%d\n" %(stage,cluster,"noPAS",x_order))




window = 60 
# canonNTcounts = numpy.zeros((window * 4 + 1,4))
# noncanonNTcounts = numpy.zeros((window * 4 + 1,4))
# noPAScounts = numpy.zeros((window * 4 + 1,4))
# canonCount = 0
# noncanonCount = 0
# noPASCount = 0
noncanonPASTable = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')
outfile.write("stage\tcluster\tPAS_type\tx_order\n")
altPAStable = []
# canon_offsets = []
# noncanon_offsets = []
for line in noncanonPASTable:
    fields = line.strip().split('\t')
    altPAStable.append(fields[0])
    

getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/L1_trimmed.fa"),outfile,altPAStable,"L1",0)
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/L2_trimmed.fa"),outfile,altPAStable,"L2",1)
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/L3_trimmed.fa"),outfile,altPAStable,"L3",2)
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/L4_trimmed.fa"),outfile,altPAStable,"L4",3)
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/YA_trimmed.fa"),outfile,altPAStable,"young adult",4)
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/GA_trimmed.fa"),outfile,altPAStable,"mature adult",5)
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/ML_trimmed.fa"),outfile,altPAStable,"male",6)
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/all_isoforms_trimmed.fa"),outfile,altPAStable,"all",7)
# getPASBreakdown(open("fastas/MG_trimmed.fa"),outfile,altPAStable,"Mangone",8)
# getPASBreakdown(open("fastas/BT_trimmed.fa"),outfile,altPAStable,"Bartel",9)
# for x in range(window * 4 + 1):
#     if sum(canonNTcounts[x,:]) != 0:
#         canonNTcounts[x,:] = 100 * canonNTcounts[x,:] / sum(canonNTcounts[x,:]) # normalize each position to a percentage
#     if sum(noncanonNTcounts[x,:]) != 0:
#         noncanonNTcounts[x,:] = 100 * noncanonNTcounts[x,:] / sum(noncanonNTcounts[x,:])
#     if sum(noPAScounts[x,:]) != 0:
#         noPAScounts[x,:] = 100 * noPAScounts[x,:] / sum(noPAScounts[x,:])
#
#
# print canonCount
# print 100.*float(canonCount) / float(canonCount+noncanonCount+noPASCount)
# print noncanonCount
# print 100.*float(noncanonCount) / float(canonCount+noncanonCount+noPASCount)
# print noPASCount
# print 100.*float(noPASCount) / float(canonCount+noncanonCount+noPASCount)
