#!/usr/bin/env python2
import sys
import matplotlib
import matplotlib.pyplot as plt
import seaborn as sns

# l1_in =  open("../../../results/utrs/beds/L1_utrs.bed")
# l2_in =  open("../../../results/utrs/beds/L2_utrs.bed")
# l3_in =  open("../../../results/utrs/beds/L3_utrs.bed")
# l4_in =  open("../../../results/utrs/beds/L4_utrs.bed")
# ya_in =  open("../../../results/utrs/beds/YA_utrs.bed")
# ga_in =  open("../../../results/utrs/beds/GA_utrs.bed")
# ml_in =  open("../../../results/utrs/beds/ML_utrs.bed")
# all_in = open("../../../results/utrs/beds/all_isoforms_utrs.bed")
#
# outfile = open("../../../results/scratch/UTRlength/utr_lengths.txt",'w')
#
# outfile.write("stage\tutr_length\tx_order\n")
# for line in l1_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("L1\t%d\t0\n" %(abs(start - end)))
# for line in l2_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("L2\t%d\t1\n" %(abs(start - end)))
#
#
# for line in l3_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("L3\t%d\t2\n" %(abs(start - end)))
#
# for line in l4_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("L4\t%d\t3\n" %(abs(start - end)))
#
# for line in ya_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("young adult\t%d\t4\n" %(abs(start - end)))
#
#
# for line in ga_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("mature adult\t%d\t5\n" %(abs(start - end)))
#
#
# for line in ml_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("male\t%d\t6\n" %(abs(start - end)))
#
# for line in all_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("all\t%d\t7\n" %(abs(start - end)))
#

### Sensitive

# l1_in =  open("../../../results/utrs/beds/L1_sensitive_utrs.bed")
# l2_in =  open("../../../results/utrs/beds/L2_sensitive_utrs.bed")
# l3_in =  open("../../../results/utrs/beds/L3_sensitive_utrs.bed")
# l4_in =  open("../../../results/utrs/beds/L4_sensitive_utrs.bed")
# ya_in =  open("../../../results/utrs/beds/YA_sensitive_utrs.bed")
# ga_in =  open("../../../results/utrs/beds/GA_sensitive_utrs.bed")
# ml_in =  open("../../../results/utrs/beds/ML_sensitive_utrs.bed")
# all_in = open("../../../results/utrs/beds/all_sensitive_isoforms_utrs.bed")
#
# outfile = open("../../../results/scratch/UTRlength/sensitive_utr_lengths.txt",'w')
#
# outfile.write("stage\tutr_length\tx_order\n")
# for line in l1_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("L1\t%d\t0\n" %(abs(start - end)))
# for line in l2_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("L2\t%d\t1\n" %(abs(start - end)))
#
#
# for line in l3_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("L3\t%d\t2\n" %(abs(start - end)))
#
# for line in l4_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("L4\t%d\t3\n" %(abs(start - end)))
#
# for line in ya_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("young adult\t%d\t4\n" %(abs(start - end)))
#
#
# for line in ga_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("mature adult\t%d\t5\n" %(abs(start - end)))
#
#
# for line in ml_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("male\t%d\t6\n" %(abs(start - end)))
#
# for line in all_in:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end = int(fields[2])
#     outfile.write("all\t%d\t7\n" %(abs(start - end)))
    
### Stringent

l1_in =  open("../../../results/utrs/beds/L1_stringent_utrs.bed")
l2_in =  open("../../../results/utrs/beds/L2_stringent_utrs.bed")
l3_in =  open("../../../results/utrs/beds/L3_stringent_utrs.bed")
l4_in =  open("../../../results/utrs/beds/L4_stringent_utrs.bed")
ya_in =  open("../../../results/utrs/beds/YA_stringent_utrs.bed")
ga_in =  open("../../../results/utrs/beds/GA_stringent_utrs.bed")
ml_in =  open("../../../results/utrs/beds/ML_stringent_utrs.bed")
all_in = open("../../../results/utrs/beds/all_stringent_isoforms_utrs.bed")

outfile = open("../../../results/scratch/UTRlength/stringent_utr_lengths.txt",'w')

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
    