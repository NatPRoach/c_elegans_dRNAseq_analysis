#!/usr/bin/env python2
import sys
import matplotlib.pyplot as plt
import seaborn as sns

# def loadspliced(read_id_to_length):
#     tmp = []
#     infile = open("../../data/L1/combined/L1_full_length.no_intron.bed")
#     for line in infile:
#         fields = line.strip().split()
#         read_id = fields[3]
#         if read_id in read_id_to_length:
#             tmp.append(read_id_to_length[read_id])
#     infile = open("../../data/L2/combined/L2_full_length.no_intron.bed")
#     for line in infile:
#         fields = line.strip().split()
#         read_id = fields[3]
#         if read_id in read_id_to_length:
#             tmp.append(read_id_to_length[read_id])
#     infile = open("../../data/L3/combined/L3_full_length.no_intron.bed")
#     for line in infile:
#         fields = line.strip().split()
#         read_id = fields[3]
#         if read_id in read_id_to_length:
#             tmp.append(read_id_to_length[read_id])
#     infile = open("../../data/L4/combined/L4_full_length.no_intron.bed")
#     for line in infile:
#         fields = line.strip().split()
#         read_id = fields[3]
#         if read_id in read_id_to_length:
#             tmp.append(read_id_to_length[read_id])
#
#     infile = open("../../data/young_adult/combined/YA_full_length.no_intron.bed")
#     for line in infile:
#         fields = line.strip().split()
#         read_id = fields[3]
#         if read_id in read_id_to_length:
#             tmp.append(read_id_to_length[read_id])
#     infile = open("../../data/adult/combined/GA_full_length.no_intron.bed")
#     for line in infile:
#         fields = line.strip().split()
#         read_id = fields[3]
#         if read_id in read_id_to_length:
#             tmp.append(read_id_to_length[read_id])
#     infile = open("../../data/male/combined/ML_full_length.no_intron.bed")
#     for line in infile:
#         fields = line.strip().split()
#         read_id = fields[3]
#         if read_id in read_id_to_length:
#             tmp.append(read_id_to_length[read_id])
#     return tmp
    
def load_sensitive_spliced(read_id_to_length):
    tmp = []
    infile = open("../../data/L1/combined/L1_sensitive_full_length.no_intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/L2/combined/L2_sensitive_full_length.no_intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/L3/combined/L3_sensitive_full_length.no_intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/L4/combined/L4_sensitive_full_length.no_intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    
    infile = open("../../data/young_adult/combined/YA_sensitive_full_length.no_intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/adult/combined/GA_sensitive_full_length.no_intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/male/combined/ML_sensitive_full_length.no_intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])  
    return tmp
    

def load_stringent_spliced(read_id_to_length):
    tmp = []
    infile = open("../../data/L1/combined/L1_stringent_full_length.no_intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/L2/combined/L2_stringent_full_length.no_intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/L3/combined/L3_stringent_full_length.no_intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/L4/combined/L4_stringent_full_length.no_intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    
    infile = open("../../data/young_adult/combined/YA_stringent_full_length.no_intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/adult/combined/GA_stringent_full_length.no_intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/male/combined/ML_stringent_full_length.no_intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])  
    return tmp

    
# def loadintron(read_id_to_length):
#     tmp = []
#     infile = open("../../data/L1/combined/L1_full_length.intron.bed")
#     for line in infile:
#         fields = line.strip().split()
#         read_id = fields[3]
#         if read_id in read_id_to_length:
#             tmp.append(read_id_to_length[read_id])
#     infile = open("../../data/L2/combined/L2_full_length.intron.bed")
#     for line in infile:
#         fields = line.strip().split()
#         read_id = fields[3]
#         if read_id in read_id_to_length:
#             tmp.append(read_id_to_length[read_id])
#     infile = open("../../data/L3/combined/L3_full_length.intron.bed")
#     for line in infile:
#         fields = line.strip().split()
#         read_id = fields[3]
#         if read_id in read_id_to_length:
#             tmp.append(read_id_to_length[read_id])
#     infile = open("../../data/L4/combined/L4_full_length.intron.bed")
#     for line in infile:
#         fields = line.strip().split()
#         read_id = fields[3]
#         if read_id in read_id_to_length:
#             tmp.append(read_id_to_length[read_id])
#
#     infile = open("../../data/young_adult/combined/YA_full_length.intron.bed")
#     for line in infile:
#         fields = line.strip().split()
#         read_id = fields[3]
#         if read_id in read_id_to_length:
#             tmp.append(read_id_to_length[read_id])
#     infile = open("../../data/adult/combined/GA_full_length.intron.bed")
#     for line in infile:
#         fields = line.strip().split()
#         read_id = fields[3]
#         if read_id in read_id_to_length:
#             tmp.append(read_id_to_length[read_id])
#     infile = open("../../data/male/combined/ML_full_length.intron.bed")
#     for line in infile:
#         fields = line.strip().split()
#         read_id = fields[3]
#         if read_id in read_id_to_length:
#             tmp.append(read_id_to_length[read_id])
#     return tmp

def load_sensitive_intron(read_id_to_length):
    tmp = []
    infile = open("../../data/L1/combined/L1_sensitive_full_length.intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/L2/combined/L2_sensitive_full_length.intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/L3/combined/L3_sensitive_full_length.intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/L4/combined/L4_sensitive_full_length.intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    
    infile = open("../../data/young_adult/combined/YA_sensitive_full_length.intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/adult/combined/GA_sensitive_full_length.intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/male/combined/ML_sensitive_full_length.intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])  
    return tmp

def load_stringent_intron(read_id_to_length):
    tmp = []
    infile = open("../../data/L1/combined/L1_stringent_full_length.intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/L2/combined/L2_stringent_full_length.intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/L3/combined/L3_stringent_full_length.intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/L4/combined/L4_stringent_full_length.intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    
    infile = open("../../data/young_adult/combined/YA_stringent_full_length.intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/adult/combined/GA_stringent_full_length.intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])
    infile = open("../../data/male/combined/ML_stringent_full_length.intron.bed")
    for line in infile:
        fields = line.strip().split()
        read_id = fields[3]
        if read_id in read_id_to_length:
            tmp.append(read_id_to_length[read_id])  
    return tmp


def loaddict():
    read_id_to_length = {}### For each read in the polyA file assign its polyA tail length in a dictionary
    tmppolyafile = open("../../data/L1/bio1/tech1/analysis/ce11_gen_L1_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L1/bio1/tech2/analysis/ce11_gen_L1_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L2/bio1/tech1/analysis/ce11_gen_L2_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L2/bio1/tech2/analysis/ce11_gen_L2_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])            

    tmppolyafile = open("../../data/L3/bio1/tech1/analysis/ce11_gen_L3_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L3/bio1/tech2/analysis/ce11_gen_L3_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L4/bio1/tech1/analysis/ce11_gen_L4_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/L4/bio1/tech2/analysis/ce11_gen_L4_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])
            
    tmppolyafile = open("../../data/young_adult/bio1/tech1/analysis/ce11_gen_young_adult_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/young_adult/bio1/tech2/analysis/ce11_gen_young_adult_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/adult/bio1/tech1/analysis/ce11_gen_adult_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/adult/bio1/tech2/analysis/ce11_gen_adult_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])
            
    tmppolyafile = open("../../data/male/bio1/tech1/analysis/ce11_gen_male_bio1_tech1_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])

    tmppolyafile = open("../../data/male/bio1/tech2/analysis/ce11_gen_male_bio1_tech2_Minimap2_k14.polya",'r')
    for i, line in enumerate(tmppolyafile):
        if i != 0:
            fields = line.strip().split()
            if fields[9] == "PASS":
                read_id_to_length[fields[0]] = float(fields[8])
    return read_id_to_length
    
read_id_to_length = loaddict()

# spliced = loadspliced(read_id_to_length)
# intron_retained = loadintron(read_id_to_length)
#
# for splice in spliced:
#     print "%0.6f\t%s" %(splice,"spliced")
# for intron in intron_retained:
#     print "%0.6f\t%s" %(intron,"intron retained")
#

# ### Sensitive
#
# spliced = load_sensitive_spliced(read_id_to_length)
# intron_retained = load_sensitive_intron(read_id_to_length)
# outfile1 = open(sys.argv[1],'w')
#
# for splice in spliced:
#     # print "%0.6f\t%s" %(splice,"spliced")
#     outfile1.write("%0.6f\t%s\n" %(splice,"spliced"))
# for intron in intron_retained:
#     # print "%0.6f\t%s" %(intron,"intron retained")
#     outfile1.write("%0.6f\t%s\n" %(intron,"intron retained"))
#
# ### Stringent
# spliced = load_stringent_spliced(read_id_to_length)
# intron_retained = load_stringent_intron(read_id_to_length)
# outfile2 = open(sys.argv[2],'w')
#
# for splice in spliced:
#     # print "%0.6f\t%s" %(splice,"spliced")
#     outfile2.write("%0.6f\t%s\n" %(splice,"spliced"))
# for intron in intron_retained:
#     # print "%0.6f\t%s" %(intron,"intron retained")
#     outfile2.write("%0.6f\t%s\n" %(intron,"intron retained"))

### Stringent
spliced = load_stringent_spliced(read_id_to_length)
intron_retained = load_stringent_intron(read_id_to_length)
outfile2 = open(sys.argv[1],'w')

for splice in spliced:
    # print "%0.6f\t%s" %(splice,"spliced")
    outfile2.write("%0.6f\t%s\n" %(splice,"spliced"))
for intron in intron_retained:
    # print "%0.6f\t%s" %(intron,"intron retained")
    outfile2.write("%0.6f\t%s\n" %(intron,"intron retained"))
