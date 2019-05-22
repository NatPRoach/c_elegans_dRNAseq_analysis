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
noncanonPASTable = open(sys.argv[1],'r')
outfile = open(sys.argv[2],'w')
outfile.write("stage\tcluster\tPAS_type\tx_order\n")
altPAStable = []
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