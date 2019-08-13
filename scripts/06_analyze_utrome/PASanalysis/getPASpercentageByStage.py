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

altPAStable = []
for line in noncanonPASTable:
    fields = line.strip().split('\t')
    altPAStable.append(fields[0])
    
outfile1 = open(sys.argv[2],'w')
outfile1.write("stage\tcluster\tPAS_type\tx_order\n")
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/L1_stringent_trimmed.fa"),outfile1,altPAStable,"L1",0)
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/L2_stringent_trimmed.fa"),outfile1,altPAStable,"L2",1)
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/L3_stringent_trimmed.fa"),outfile1,altPAStable,"L3",2)
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/L4_stringent_trimmed.fa"),outfile1,altPAStable,"L4",3)
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/YA_stringent_trimmed.fa"),outfile1,altPAStable,"young adult",4)
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/GA_stringent_trimmed.fa"),outfile1,altPAStable,"mature adult",5)
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/ML_stringent_trimmed.fa"),outfile1,altPAStable,"male",6)
getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/all_stringent_isoforms_trimmed.fa"),outfile1,altPAStable,"all",7)

# outfile2 = open(sys.argv[3],'w')
# outfile2.write("stage\tcluster\tPAS_type\tx_order\n")
# getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/L1_sensitive_trimmed.fa"),outfile2,altPAStable,"L1",0)
# getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/L2_sensitive_trimmed.fa"),outfile2,altPAStable,"L2",1)
# getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/L3_sensitive_trimmed.fa"),outfile2,altPAStable,"L3",2)
# getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/L4_sensitive_trimmed.fa"),outfile2,altPAStable,"L4",3)
# getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/YA_sensitive_trimmed.fa"),outfile2,altPAStable,"young adult",4)
# getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/GA_sensitive_trimmed.fa"),outfile2,altPAStable,"mature adult",5)
# getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/ML_sensitive_trimmed.fa"),outfile2,altPAStable,"male",6)
# getPASBreakdown(open("../../../results/scratch/PASanalysis/fastas/all_sensitive_isoforms_trimmed.fa"),outfile2,altPAStable,"all",7)