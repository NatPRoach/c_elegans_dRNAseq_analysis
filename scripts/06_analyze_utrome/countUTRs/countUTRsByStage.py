#!/usr/bin/env python2
import sets
l1_mangone_in  = open("../../../references/utrs/l1_mangone_utrs.bed")
l2_mangone_in  = open("../../../references/utrs/l2_mangone_utrs.bed")
l3_mangone_in  = open("../../../references/utrs/l3_mangone_utrs.bed")
l4_mangone_in  = open("../../../references/utrs/l4_mangone_utrs.bed")
ya_mangone_in  = open("../../../references/utrs/ya_mangone_utrs.bed")
ga_mangone_in  = open("../../../references/utrs/ya_mangone_utrs.bed")
ml_mangone_in  = open("../../../references/utrs/ml_mangone_utrs.bed")
all_mangone_in = open("../../../references/utrs/mangone_utrs.bed")

l1_mangone_uniq = sets.Set()
l2_mangone_uniq = sets.Set()
l3_mangone_uniq = sets.Set()
l4_mangone_uniq = sets.Set()
ya_mangone_uniq = sets.Set()
ga_mangone_uniq = sets.Set()
ml_mangone_uniq = sets.Set()
all_mangone_uniq = sets.Set()

for line in l1_mangone_in:
    fields = line.strip().split()
    l1_mangone_uniq.add((fields[0],fields[1],fields[2],fields[5]))
for line in l2_mangone_in:
    fields = line.strip().split()
    l2_mangone_uniq.add((fields[0],fields[1],fields[2],fields[5]))
for line in l3_mangone_in:
    fields = line.strip().split()
    l3_mangone_uniq.add((fields[0],fields[1],fields[2],fields[5]))
for line in l4_mangone_in:
    fields = line.strip().split()
    l4_mangone_uniq.add((fields[0],fields[1],fields[2],fields[5]))
for line in ya_mangone_in:
    fields = line.strip().split()
    ya_mangone_uniq.add((fields[0],fields[1],fields[2],fields[5]))
for line in ga_mangone_in:
    fields = line.strip().split()
    ga_mangone_uniq.add((fields[0],fields[1],fields[2],fields[5]))
for line in ml_mangone_in:
    fields = line.strip().split()
    ml_mangone_uniq.add((fields[0],fields[1],fields[2],fields[5]))
for line in all_mangone_in:
    fields = line.strip().split()
    all_mangone_uniq.add((fields[0],fields[1],fields[2],fields[5]))


all_jan_in = open("../../../references/utrs/jan_utrs.bed")

all_jan_uniq = sets.Set()
for line in all_jan_in:
    fields = line.strip().split()
    all_jan_uniq.add((fields[0],fields[1],fields[2],fields[5]))

# l1_in  = open("../../../results/utrs/beds/L1_utrs.bed")
# l2_in  = open("../../../results/utrs/beds/L2_utrs.bed")
# l3_in  = open("../../../results/utrs/beds/L3_utrs.bed")
# l4_in  = open("../../../results/utrs/beds/L4_utrs.bed")
# ya_in  = open("../../../results/utrs/beds/YA_utrs.bed")
# ga_in  = open("../../../results/utrs/beds/GA_utrs.bed")
# ml_in  = open("../../../results/utrs/beds/ML_utrs.bed")
# all_in = open("../../../results/utrs/beds/all_isoforms_utrs.bed")
#
#
# l1_uniq = sets.Set()
# l2_uniq = sets.Set()
# l3_uniq = sets.Set()
# l4_uniq = sets.Set()
# ya_uniq = sets.Set()
# ga_uniq = sets.Set()
# ml_uniq = sets.Set()
# all_uniq = sets.Set()
# for line in l1_in:
#     fields = line.strip().split()
#     l1_uniq.add((fields[0],fields[1],fields[2],fields[5]))
# for line in l2_in:
#     fields = line.strip().split()
#     l2_uniq.add((fields[0],fields[1],fields[2],fields[5]))
# for line in l3_in:
#     fields = line.strip().split()
#     l3_uniq.add((fields[0],fields[1],fields[2],fields[5]))
# for line in l4_in:
#     fields = line.strip().split()
#     l4_uniq.add((fields[0],fields[1],fields[2],fields[5]))
# for line in ya_in:
#     fields = line.strip().split()
#     ya_uniq.add((fields[0],fields[1],fields[2],fields[5]))
# for line in ga_in:
#     fields = line.strip().split()
#     ga_uniq.add((fields[0],fields[1],fields[2],fields[5]))
# for line in ml_in:
#     fields = line.strip().split()
#     ml_uniq.add((fields[0],fields[1],fields[2],fields[5]))
# for line in all_in:
#     fields = line.strip().split()
#     all_uniq.add((fields[0],fields[1],fields[2],fields[5]))
#
#
#
# outfile = open("../../../results/scratch/countUTRs/utr_counts.txt",'w')
#
#
# outfile.write("stage\tdataset\tcounts\tx_order\n")
#
# outfile.write("L1\tOur UTRs\t%d\t0\n" %(len(l1_uniq)))
# outfile.write("L2\tOur UTRs\t%d\t1\n" %(len(l2_uniq)))
# outfile.write("L3\tOur UTRs\t%d\t2\n" %(len(l3_uniq)))
# outfile.write("L4\tOur UTRs\t%d\t3\n" %(len(l4_uniq)))
# outfile.write("young adult\tOur UTRs\t%d\t4\n" %(len(ya_uniq)))
# outfile.write("mature adult\tOur UTRs\t%d\t5\n" %(len(ga_uniq)))
# outfile.write("male\tOur UTRs\t%d\t6\n" %(len(ml_uniq)))
# outfile.write("all\tOur UTRs\t%d\t7\n" %(len(all_uniq)))
#
# outfile.write("L1\tMangone UTRs\t%d\t0\n" %(len(l1_mangone_uniq)))
# outfile.write("L2\tMangone UTRs\t%d\t1\n" %(len(l2_mangone_uniq)))
# outfile.write("L3\tMangone UTRs\t%d\t2\n" %(len(l3_mangone_uniq)))
# outfile.write("L4\tMangone UTRs\t%d\t3\n" %(len(l4_mangone_uniq)))
# outfile.write("young adult\tMangone UTRs\t%d\t4\n" %(len(ya_mangone_uniq)))
# outfile.write("mature adult\tMangone UTRs\t%d\t5\n" %(len(ga_mangone_uniq)))
# outfile.write("male\tMangone UTRs\t%d\t6\n" %(len(ml_mangone_uniq)))
# outfile.write("all\tMangone UTRs\t%d\t7\n" %(len(all_mangone_uniq)))
#
# outfile.write("all\tJan UTRs\t%d\t7\n" %(len(all_jan_uniq)))

# ### Sensitive
#
# l1_in  = open("../../../results/utrs/beds/L1_sensitive_utrs.bed")
# l2_in  = open("../../../results/utrs/beds/L2_sensitive_utrs.bed")
# l3_in  = open("../../../results/utrs/beds/L3_sensitive_utrs.bed")
# l4_in  = open("../../../results/utrs/beds/L4_sensitive_utrs.bed")
# ya_in  = open("../../../results/utrs/beds/YA_sensitive_utrs.bed")
# ga_in  = open("../../../results/utrs/beds/GA_sensitive_utrs.bed")
# ml_in  = open("../../../results/utrs/beds/ML_sensitive_utrs.bed")
# all_in = open("../../../results/utrs/beds/all_sensitive_isoforms_utrs.bed")
#
#
# l1_uniq = sets.Set()
# l2_uniq = sets.Set()
# l3_uniq = sets.Set()
# l4_uniq = sets.Set()
# ya_uniq = sets.Set()
# ga_uniq = sets.Set()
# ml_uniq = sets.Set()
# all_uniq = sets.Set()
# for line in l1_in:
#     fields = line.strip().split()
#     l1_uniq.add((fields[0],fields[1],fields[2],fields[5]))
# for line in l2_in:
#     fields = line.strip().split()
#     l2_uniq.add((fields[0],fields[1],fields[2],fields[5]))
# for line in l3_in:
#     fields = line.strip().split()
#     l3_uniq.add((fields[0],fields[1],fields[2],fields[5]))
# for line in l4_in:
#     fields = line.strip().split()
#     l4_uniq.add((fields[0],fields[1],fields[2],fields[5]))
# for line in ya_in:
#     fields = line.strip().split()
#     ya_uniq.add((fields[0],fields[1],fields[2],fields[5]))
# for line in ga_in:
#     fields = line.strip().split()
#     ga_uniq.add((fields[0],fields[1],fields[2],fields[5]))
# for line in ml_in:
#     fields = line.strip().split()
#     ml_uniq.add((fields[0],fields[1],fields[2],fields[5]))
# for line in all_in:
#     fields = line.strip().split()
#     all_uniq.add((fields[0],fields[1],fields[2],fields[5]))
#
#
#
# outfile = open("../../../results/scratch/countUTRs/sensitive_utr_counts.txt",'w')
#
#
# outfile.write("stage\tdataset\tcounts\tx_order\n")
#
# outfile.write("L1\tOur UTRs\t%d\t0\n" %(len(l1_uniq)))
# outfile.write("L2\tOur UTRs\t%d\t1\n" %(len(l2_uniq)))
# outfile.write("L3\tOur UTRs\t%d\t2\n" %(len(l3_uniq)))
# outfile.write("L4\tOur UTRs\t%d\t3\n" %(len(l4_uniq)))
# outfile.write("young adult\tOur UTRs\t%d\t4\n" %(len(ya_uniq)))
# outfile.write("mature adult\tOur UTRs\t%d\t5\n" %(len(ga_uniq)))
# outfile.write("male\tOur UTRs\t%d\t6\n" %(len(ml_uniq)))
# outfile.write("all\tOur UTRs\t%d\t7\n" %(len(all_uniq)))
#
# outfile.write("L1\tMangone UTRs\t%d\t0\n" %(len(l1_mangone_uniq)))
# outfile.write("L2\tMangone UTRs\t%d\t1\n" %(len(l2_mangone_uniq)))
# outfile.write("L3\tMangone UTRs\t%d\t2\n" %(len(l3_mangone_uniq)))
# outfile.write("L4\tMangone UTRs\t%d\t3\n" %(len(l4_mangone_uniq)))
# outfile.write("young adult\tMangone UTRs\t%d\t4\n" %(len(ya_mangone_uniq)))
# outfile.write("mature adult\tMangone UTRs\t%d\t5\n" %(len(ga_mangone_uniq)))
# outfile.write("male\tMangone UTRs\t%d\t6\n" %(len(ml_mangone_uniq)))
# outfile.write("all\tMangone UTRs\t%d\t7\n" %(len(all_mangone_uniq)))
#
# outfile.write("all\tJan UTRs\t%d\t7\n" %(len(all_jan_uniq)))


### Stringent

l1_in  = open("../../../results/utrs/beds/L1_stringent_utrs.bed")
l2_in  = open("../../../results/utrs/beds/L2_stringent_utrs.bed")
l3_in  = open("../../../results/utrs/beds/L3_stringent_utrs.bed")
l4_in  = open("../../../results/utrs/beds/L4_stringent_utrs.bed")
ya_in  = open("../../../results/utrs/beds/YA_stringent_utrs.bed")
ga_in  = open("../../../results/utrs/beds/GA_stringent_utrs.bed")
ml_in  = open("../../../results/utrs/beds/ML_stringent_utrs.bed")
all_in = open("../../../results/utrs/beds/all_stringent_isoforms_utrs.bed")


l1_uniq = sets.Set()
l2_uniq = sets.Set()
l3_uniq = sets.Set()
l4_uniq = sets.Set()
ya_uniq = sets.Set()
ga_uniq = sets.Set()
ml_uniq = sets.Set()
all_uniq = sets.Set()
for line in l1_in:
    fields = line.strip().split()
    l1_uniq.add((fields[0],fields[1],fields[2],fields[5]))
for line in l2_in:
    fields = line.strip().split()
    l2_uniq.add((fields[0],fields[1],fields[2],fields[5]))
for line in l3_in:
    fields = line.strip().split()
    l3_uniq.add((fields[0],fields[1],fields[2],fields[5]))
for line in l4_in:
    fields = line.strip().split()
    l4_uniq.add((fields[0],fields[1],fields[2],fields[5]))
for line in ya_in:
    fields = line.strip().split()
    ya_uniq.add((fields[0],fields[1],fields[2],fields[5]))
for line in ga_in:
    fields = line.strip().split()
    ga_uniq.add((fields[0],fields[1],fields[2],fields[5]))
for line in ml_in:
    fields = line.strip().split()
    ml_uniq.add((fields[0],fields[1],fields[2],fields[5]))
for line in all_in:
    fields = line.strip().split()
    all_uniq.add((fields[0],fields[1],fields[2],fields[5]))



outfile = open("../../../results/scratch/countUTRs/stringent_utr_counts.txt",'w')


outfile.write("stage\tdataset\tcounts\tx_order\n")

outfile.write("L1\tOur UTRs\t%d\t0\n" %(len(l1_uniq)))
outfile.write("L2\tOur UTRs\t%d\t1\n" %(len(l2_uniq)))
outfile.write("L3\tOur UTRs\t%d\t2\n" %(len(l3_uniq)))
outfile.write("L4\tOur UTRs\t%d\t3\n" %(len(l4_uniq)))
outfile.write("young adult\tOur UTRs\t%d\t4\n" %(len(ya_uniq)))
outfile.write("mature adult\tOur UTRs\t%d\t5\n" %(len(ga_uniq)))
outfile.write("male\tOur UTRs\t%d\t6\n" %(len(ml_uniq)))
outfile.write("all\tOur UTRs\t%d\t7\n" %(len(all_uniq)))

outfile.write("L1\tMangone UTRs\t%d\t0\n" %(len(l1_mangone_uniq)))
outfile.write("L2\tMangone UTRs\t%d\t1\n" %(len(l2_mangone_uniq)))
outfile.write("L3\tMangone UTRs\t%d\t2\n" %(len(l3_mangone_uniq)))
outfile.write("L4\tMangone UTRs\t%d\t3\n" %(len(l4_mangone_uniq)))
outfile.write("young adult\tMangone UTRs\t%d\t4\n" %(len(ya_mangone_uniq)))
outfile.write("mature adult\tMangone UTRs\t%d\t5\n" %(len(ga_mangone_uniq)))
outfile.write("male\tMangone UTRs\t%d\t6\n" %(len(ml_mangone_uniq)))
outfile.write("all\tMangone UTRs\t%d\t7\n" %(len(all_mangone_uniq)))

outfile.write("all\tJan UTRs\t%d\t7\n" %(len(all_jan_uniq)))