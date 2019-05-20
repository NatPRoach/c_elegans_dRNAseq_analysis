#!/usr/bin/env python2
import sets

def parseStage(info):
    gene_id,source,library = info.split(',')
    foo = source.split('#')
    if len(foo) != 2:
        return None
    s1, s2 = foo
    if s2 == "RNA":
        return s1
    elif s2 == "Capture":
        if s1[:4] == "U454":
            return None
        else:
            return s1.split("_")[0]
    elif s2 == "Clone" and library != "NULL":
        if library == "ycL1":
            return "L1"
        elif library == "ycL2":
            return "L2"
        elif library == "ycL3":
            return "L3"
        elif library == "ycL4":
            return "L4"
        else:
            return None
    else:
        return None


infile = open("../../references/utrs/mangone_utrs.bed")
l1_out = open("../../references/utrs/l1_mangone_utrs.bed",'w')
l2_out = open("../../references/utrs/l2_mangone_utrs.bed",'w')
l3_out = open("../../references/utrs/l3_mangone_utrs.bed",'w')
l4_out = open("../../references/utrs/l4_mangone_utrs.bed",'w')
ya_out = open("../../references/utrs/ya_mangone_utrs.bed",'w')
# ga_out = open("ga_mangone_utrs.bed",'w')
ml_out = open("../../references/utrs/ml_mangone_utrs.bed",'w')

l1_uniq = sets.Set()
l2_uniq = sets.Set()
l3_uniq = sets.Set()
l4_uniq = sets.Set()
ya_uniq = sets.Set()
# ga_uniq = sets.Set()
ml_uniq = sets.Set()

for line in infile:
    fields = line.strip().split()
    info = fields[3]
    stage = parseStage(info)
    if stage is None:
        continue
    if stage == "L1":
        l1_out.write(line)
        l1_uniq.add((fields[0],fields[1],fields[2]))
    elif stage == "L2":
        l2_out.write(line)
        l2_uniq.add((fields[0],fields[1],fields[2]))
    elif stage == "L3":
        l3_out.write(line)
        l3_uniq.add((fields[0],fields[1],fields[2]))
    elif stage == "L4":
        l4_out.write(line)
        l4_uniq.add((fields[0],fields[1],fields[2]))
    elif stage == "YA":
        ya_out.write(line)
        ya_uniq.add((fields[0],fields[1],fields[2]))
    elif stage == "male":
        ml_out.write(line)
        ml_uniq.add((fields[0],fields[1],fields[2]))
    elif stage == "herm":
        ya_out.write(line)
        ya_uniq.add((fields[0],fields[1],fields[2]))
    elif stage == "emb":
        continue
print len(l1_uniq)
print len(l2_uniq)
print len(l3_uniq)
print len(l4_uniq)
print len(ya_uniq)
# print len(ga_uniq)
print len(ml_uniq)