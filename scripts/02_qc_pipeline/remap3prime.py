#!/usr/bin/env python2

import sys
import numpy as np
import pysam
import subprocess
import os
import tempfile
def s(x,y,a,b):
    if x == y:
        return a
    else:
        return -b
def align(x,y,a=1,b=2,o=2,e=1,anchor=None):
    #x - first seq (query for the purposes of defining insertion / deletion)
    #y - second seq (reference for the purposes of defining insertion / deletion)
    #a - match score
    #b - mismatch penalty
    #o - gap open penalty
    #e - gap extend penalty
    inf = np.inf
    x = x.upper()
    y = y.upper()
    v = np.zeros([len(x) + 1, len(y) + 1])
    m = np.zeros([len(x) + 1, len(y) + 1])
    x_gap = np.zeros([len(x) + 1, len(y) + 1])
    y_gap = np.zeros([len(x) + 1, len(y) + 1])
    walkback_matrix = np.zeros([len(x) + 1, len(y) + 1])
    
    ### Initialize matrix values
    if anchor is None: #Semi-global alignment, allow gaps at the begining and end of y
        ## Both end gaps are free. (for y only)
        for i in range(1,len(x)+1): 
            v[i,0] = -o - i*e
            #v[i,0] = 0. #Redundant, but do it for now
            #x_gap[i,0] = -o - i*e
            y_gap[i,0] = -o - i*e
            #y_gap[i,0] = -inf
            #x_gap[i,0] = -inf
        for j in range(1,len(y)+1): 
            #v[0,j] = -o - j*e
            x_gap[0,j] = -o - j*e
            #y_gap[0,j] = -o - j*e
            v[0,j] = 0. #Redundant, but do it for now
            #x_gap[0,j] = -inf
            #y_gap[0,j] = -inf
    elif anchor == "both": ## Working as intended
        ## Both end gaps are penalized
        for i in range(1,len(x)+1): 
            v[i,0] = -o - i*e
            x_gap[i,0] = -o - i*e
            y_gap[i,0] = -inf
        for j in range(1,len(y)+1):
            v[0,j] = -o - j*e
            y_gap[0,j] = -o - j*e
            x_gap[0,j] = -inf
    elif anchor == "left" or anchor == "l" or anchor == "L":
        ## Gaps on the left side are penalized
        for i in range(1,len(x)+1): 
            v[i,0] = -o - i*e
            #v[i,0] = 0. #Redundant, but do it for now
            x_gap[i,0] = -o - i*e
            #y_gap[i,0] = -o - i*e
            y_gap[i,0] = -inf
            #x_gap[i,0] = -inf
        for j in range(1,len(y)+1): 
            v[0,j] = -o - j*e
            #x_gap[0,j] = -o - j*e
            y_gap[0,j] = -o - j*e
            #v[0,j] = 0. #Redundant, but do it for now
            x_gap[0,j] = -inf
            #y_gap[0,j] = -inf
    elif anchor == "right" or anchor == "r" or anchor == "R":
        ## Gaps on the right side are penalized
        for i in range(1,len(x)+1): 
            v[i,0] = -o - i*e
            #v[i,0] = 0. #Redundant, but do it for now
            x_gap[i,0] = -o - i*e
            #y_gap[i,0] = -o - i*e
            y_gap[i,0] = -inf
            #x_gap[i,0] = -inf
        for j in range(1,len(y)+1): 
            #v[0,j] = -o - j*e
            #x_gap[0,j] = -o - j*e
            y_gap[0,j] = -o - j*e
            v[0,j] = 0. #Redundant, but do it for now
            x_gap[0,j] = -inf
            #y_gap[0,j] = -inf
    ## Fill in the score matrices:
    for i in range(1,len(x)+1):
        for j in range(1,len(y)+1):
            m[i,j] = v[i-1,j-1] + s(x[i-1],y[j-1],a,b)
            x_gap[i,j] = max( x_gap[i-1,j], v[i-1,j] - o) - e
            y_gap[i,j] = max( y_gap[i,j-1], v[i,j-1] - o) - e
            v[i,j] = max([m[i,j],x_gap[i,j],y_gap[i,j]])
            possible_origins = 0
            if v[i,j] == m[i,j]:
                possible_origins += 1
            if v[i,j] == x_gap[i,j]:
                possible_origins += 2
            if v[i,j] == y_gap[i,j]:
                possible_origins += 4
            walkback_matrix[i,j] = possible_origins
    ## Walkback
    if anchor is None or anchor == "left":
        i = len(x)
        #print i
        optimal_score = max(v[i,:])
        #print optimal_score
        for k in range(len(y)+1):
            if v[i,k] == optimal_score:
                j = k
                break
        #print "j"
        #print j
        walkback = [ 2 for x in range(len(y)- k)] 
    elif anchor == "both" or anchor == "right":
        i = len(x)
        j = len(y)
        walkback = []
    
    while i != 0 and j != 0:
        direction = walkback_matrix[i,j]
        if direction == 1:
            i = i-1
            j = j-1
            walkback.append(0)
        elif direction == 2:
            i= i-1
            walkback.append(1)
        elif direction == 4:
            j=j-1
            walkback.append(2)
        elif direction == 3: #Ambiguous alignment
            i = i-1 #Arbitrarily decide to use match
            j = j-1
            walkback.append(0)
        elif direction == 5: #Ambiguous alignment
            i = i-1 #Arbitrarily decide to use match
            j = j-1
            walkback.append(0)
        elif direction == 6: #Ambiguous alignment
            j=j-1 #Arbitrarily decide to use deletion
            walkback.append(2)
        elif direction == 7: #Ambiguous alignment
            i = i-1 #Arbitrarily decide to use match
            j = j-1
            walkback.append(0)
    if i > 0:
        walkback =  walkback + [1 for x in range(i)]
    elif j > 0:
        walkback =  walkback + [2 for x in range(j)]
    walkback.reverse()
    cigar = convertWalkbackToCIGAR(walkback,anchor= anchor)
    return cigar

def convertWalkbackToCIGAR(w,anchor=None):
    cigar = []
    flag = False
    last_op = None
    counter = 1
    if anchor is None or anchor == "both":
        for i, c in enumerate(w):
            if not flag:
                if c == 1:
                    flag = True
                else:
                    continue
            if c == last_op:
                counter += 1
            else:
                if last_op is not None:
                    cigar.append( (last_op,counter))
                last_op = c
                counter = 1
        cigar.append((last_op,counter))
        if cigar[-1][0] == 2:
            cigar = cigar[:-1]
    elif anchor == "left":
        for c in w:
            if c == last_op:
                counter +=1
            else:
                if last_op is not None:
                    cigar.append((last_op,counter))
                last_op = c
                counter = 1
        cigar.append((last_op,counter))
        if cigar[-1][0] == 2:
            cigar = cigar[:-1]
        if cigar[-1][0] == 1:
            cigar[-1] = (4, cigar[-1][1])
    elif anchor == "right":
        for i, c in enumerate(w):
            if not flag:
                if c == 0:
                    flag = True
                else:
                    continue
            if c == last_op:
                counter += 1
            else:
                if last_op is not None:
                    cigar.append( (last_op,counter))
                last_op = c
                counter = 1
        cigar.append((last_op,counter))
        if cigar[0][0] == 1:
            cigar[0] = (4,cigar[0][1])
    return cigar

def calcRefStart(old_ref_start, cigar):
    extend = 0 
    for op,length in cigar:
        if op == 4:
            continue
        elif op == 1:
            continue
        elif op == 2:
            extend += length
        elif op == 0:
            extend += length
    return old_ref_start - extend

def trimSoftClip(seq,quals,threshold=12.,anchor="left"):
    if anchor == "left":
        flag = False
        for x in reversed(range(len(seq))):
            if seq[x] != 'A':
                flag = True
                break
        if flag: #Seq has non A char
            trim = seq[:x+1]
            soft_clip_amount = len(seq) - x - 1
            
            qual_block = []
            a_block_start = None
            flag2 = False
            for x in range(len(seq)):
                if seq[x] == "A":
                    qual_block.append(quals[x])
                    if not flag2:
                        a_block_start = x
                        flag2 = True
                else:
                    flag2 = False
                    if len(qual_block) >= 3 and np.mean(qual_block) > threshold:
                        break
                    qual_block = []
            if len(qual_block) >= 3:
                if np.mean(qual_block) > threshold:
                    trim2 = seq[:a_block_start]
                    soft_clip_amount2 = len(seq) - a_block_start
                    if soft_clip_amount < soft_clip_amount2:
                        trim = trim2
                        soft_clip_amount = soft_clip_amount2
        else:#Seq is all As
            trim = ""
            soft_clip_amount = len(seq)
    elif anchor == "right":
        flag = False
        for x in range(len(seq)):
            if seq[x] != 'T':
                flag = True
                break
        if flag: #Seq has non A char
            trim = seq[x:]
            soft_clip_amount = x
            
            qual_block = []
            a_block_start = None
            flag2 = False
            for x in reversed(range(len(seq))):
                if seq[x] == "T":
                    qual_block.append(quals[x])
                    if not flag2:
                        a_block_start = x
                        flag2 = True
                else:
                    flag2 = False
                    if len(qual_block) >= 3 and np.mean(qual_block) > threshold:
                        break
                    qual_block = []
            if len(qual_block) >= 3:
                if np.mean(qual_block) > threshold:
                    trim2 = seq[a_block_start + 1:]
                    soft_clip_amount2 = a_block_start  + 1
                    if soft_clip_amount < soft_clip_amount2:
                        trim = trim2
                        soft_clip_amount = soft_clip_amount2
            
        else:#Seq is all As
            trim = ""
            soft_clip_amount = len(seq)
    else:
        print "Error in trim poly(A)"
    return soft_clip_amount, trim

def trimPolyA(seq,anchor="left"):
    if anchor == "left":
        flag = False
        for x in reversed(range(len(seq))):
            if seq[x] != 'A':
                flag = True
                break
        if flag: #Seq has non A char
            trim = seq[:x+1]
            soft_clip_amount = len(seq) - x - 1
        else:#Seq is all As
            trim = ""
            soft_clip_amount = len(seq)
    elif anchor == "right":
        flag = False
        for x in range(len(seq)):
            if seq[x] != 'T':
                flag = True
                break
        if flag: #Seq has non A char
            trim = seq[x:]
            soft_clip_amount = x
        else:#Seq is all As
            trim = ""
            soft_clip_amount = len(seq)
    else:
        print "Error in trim poly(A)"
    return soft_clip_amount, trim

def formatCIGAR(cigar_list):
    out_cigar = []
    for op,length in cigar_list:
        out_cigar.append("%d%s" %(length,op))
    return "".join(out_cigar)

def testCIGAR(cigar_list,seq_len):
    total = 0
    for op,length in cigar_list:
        if op == 'M':
            total += length
        elif op == 'I':
            total += length
        elif op == 'S':
            total += length
        elif op == '=':
            total += length
        elif op == 'X':
            total += length
    return total == seq_len

infile = pysam.AlignmentFile(sys.argv[1],'rb')
outfile = pysam.AlignmentFile(sys.argv[2],'wb',template=infile)
exclusion_file = open(sys.argv[3],'w')
ref_in = open(sys.argv[4])

window = 300

ref = {}
for i,line in enumerate(ref_in):
    if line[0] == '>':
        if i != 0:
            ref[chrom] = seq
        seq = ""
        chrom = line.strip("\n>")
    else:
        seq = seq + line.strip()


for read in infile.fetch():
    seq = read.query_sequence
    qual = read.query_qualities
    #print qual
    chrom = read.reference_name
    if chrom in ref:
        if read.is_reverse: #reverse strand
            cigar_tuples = read.cigartuples
            if cigar_tuples[0][0] == 4: #softclipped?
                soft_clip_len = cigar_tuples[0][1]
                if soft_clip_len >= 10:
                    exclusion_file.write("%s\n" %(read.query_name))
                side = "right"
                soft_bases = seq[:soft_clip_len]
                soft_quals = qual[:soft_clip_len]
                ref_bases = ref[chrom][(read.reference_start-window):read.reference_start].upper()
                soft_clip_cigar = align(soft_bases,ref_bases,anchor="right")
                dup = soft_clip_cigar + cigar_tuples[1:]
                read.cigartuples = dup
                read.reference_start = calcRefStart(read.reference_start,soft_clip_cigar)
        else: #forward strand
            cigar_tuples = read.cigartuples
            if cigar_tuples[-1][0] == 4: #softclipped?
                soft_clip_len = cigar_tuples[-1][1]
                if soft_clip_len >= 10:
                    exclusion_file.write("%s\n" %(read.query_name))
                side = "left"
                soft_bases = seq[-soft_clip_len:]
                soft_quals = qual[-soft_clip_len:]
                ref_bases = ref[chrom][(read.reference_end):read.reference_end + window].upper()
                soft_clip_cigar = align(soft_bases,ref_bases,anchor="left")
                dup = cigar_tuples[:-1] + soft_clip_cigar
                read.cigartuples = dup
        
    outfile.write(read)
