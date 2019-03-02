#!/usr/bin/env python2
import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

l1_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/L1_utrs.bed")
l2_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/L2_utrs.bed")
l3_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/L3_utrs.bed")
l4_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/L4_utrs.bed")
ya_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/YA_utrs.bed")
ga_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/GA_utrs.bed")
ml_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/ML_utrs.bed")
all_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/all_isoforms_utrs.bed")

outfile = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/analyzeUTRome/UTRlength/utr_lengths.txt",'w')

outfile.write("stage\tutr_length\tx_order\n")
for line in l1_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile.write("L1\t%d\t0\n" %(abs(start - end)))    
for line in l2_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile.write("L2\t%d\t1\n" %(abs(start - end)))    

    
for line in l3_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile.write("L3\t%d\t2\n" %(abs(start - end)))
    
for line in l4_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile.write("L4\t%d\t3\n" %(abs(start - end)))    
    
for line in ya_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile.write("young adult\t%d\t4\n" %(abs(start - end)))    

    
for line in ga_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile.write("mature adult\t%d\t5\n" %(abs(start - end)))    

    
for line in ml_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile.write("male\t%d\t6\n" %(abs(start - end)))
    
for line in all_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile.write("all\t%d\t7\n" %(abs(start - end)))


l1_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/alls/all_L1_utrs.bed")
l2_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/alls/all_L2_utrs.bed")
l3_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/alls/all_L3_utrs.bed")
l4_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/alls/all_L4_utrs.bed")
ya_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/alls/all_YA_utrs.bed")
ga_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/alls/all_GA_utrs.bed")
ml_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/alls/all_ML_utrs.bed")
all_in = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/results/realigned_utrs/beds/alls/all_utrs.bed")

outfile2 = open("/Users/nproach/Documents/NPR_Notebook/04_Define_Isoforms/analyzeUTRome/UTRlength/all_utr_lengths.txt",'w')

outfile2.write("stage\tutr_length\tx_order\n")
for line in l1_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile2.write("L1\t%d\t0\n" %(abs(start - end)))    
for line in l2_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile2.write("L2\t%d\t1\n" %(abs(start - end)))    

    
for line in l3_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile2.write("L3\t%d\t2\n" %(abs(start - end)))
    
for line in l4_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile2.write("L4\t%d\t3\n" %(abs(start - end)))    
    
for line in ya_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile2.write("young adult\t%d\t4\n" %(abs(start - end)))    

    
for line in ga_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile2.write("mature adult\t%d\t5\n" %(abs(start - end)))    

    
for line in ml_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile2.write("male\t%d\t6\n" %(abs(start - end)))
    
for line in all_in:
    fields = line.strip().split()
    start = int(fields[1])
    end = int(fields[2])
    outfile2.write("all\t%d\t7\n" %(abs(start - end)))
