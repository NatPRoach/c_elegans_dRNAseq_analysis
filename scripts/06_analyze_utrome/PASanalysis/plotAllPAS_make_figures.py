#!/usr/bin/env python2

import sys
import numpy
import matplotlib
import matplotlib.pyplot as plt


def plotAllPAS(fastaFile,noncanonPASTable,pasAssignmentFile,outfiles_prefix):
    font = {"family":"sans-serif",
            "sans-serif":["Helvetica"],
            "size":8}
    matplotlib.rc('font',**font)
    
    window = 60 
    canonNTcounts = numpy.zeros((window * 4 + 1,4))
    noncanonNTcounts = numpy.zeros((window * 4 + 1,4))
    noPAScounts = numpy.zeros((window * 4 + 1,4))
    canonCount = 0
    noncanonCount = 0
    noPASCount = 0
#     fastaFile = open(sys.argv[1],'r')
#     noncanonPASTable = open(sys.argv[2],'r')
#     pasAssignmentFile = open(sys.argv[3],'w')

    altPAStable = []
    canon_offsets = []
    noncanon_offsets = []
    for line in noncanonPASTable:
        fields = line.strip().split('\t')
        altPAStable.append(fields[0])
    last_line = None 
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
                offset = window * 2  - 19 - x
                canon_offsets.append(x - window) ## TODO: Make sure this math is right
                canonFlag = True
                break
        if canonFlag:
            canonCount += 1
            for x in range(0,len(seq)):
                if x + offset >= 0 and x + offset < 4*window+1:
                    nt = seq[x]
                    if nt == "A":
                        canonNTcounts[x+offset,0] += 1
                    elif nt == "C":
                        canonNTcounts[x+offset,1] += 1
                    elif nt == "G":
                        canonNTcounts[x+offset,2] += 1
                    elif nt == "T":
                        canonNTcounts[x+offset,3] += 1
                    else:
                        print 'something went wrong'
            pasAssignmentFile.write("%s\t%s\n" %(cluster,"AATAAA"))
        else:
            noncanonFlag = False
            for altPAS in altPAStable:
                for x in reversed(range(20,window-5)):
                    subseq = seq[x:x+6]
                    if subseq == altPAS:
                        offset = window * 2 - 19 - x
                        noncanon_offsets.append(x - window ) ## TODO: Make sure this math is right
                        noncanonFlag = True
                        #print altPAS
                        break
                if noncanonFlag:
                    break
            if noncanonFlag:
                noncanonCount += 1
                for x in range(0,len(seq)):
                    if x + offset >= 0 and x + offset < window * 4 + 1:
                        nt = seq[x]
                        if nt == "A":
                            noncanonNTcounts[x+offset,0] += 1
                        elif nt == "C":
                            noncanonNTcounts[x+offset,1] += 1
                        elif nt == "G":
                            noncanonNTcounts[x+offset,2] += 1
                        elif nt == "T":
                            noncanonNTcounts[x+offset,3] += 1
                        else:
                            print 'something went wrong'
                pasAssignmentFile.write("%s\t%s\n" %(cluster,altPAS))
            else:
                #print last_line
                #print line.strip()
                noPASCount += 1
                offset = window
                for x in range(0,len(seq)):
                    if x + offset >= 0 and x + offset < window * 4 + 1:
                        nt = seq[x]
                        if nt == "A":
                            noPAScounts[x+offset,0] += 1
                        elif nt == "C":
                            noPAScounts[x+offset,1] += 1
                        elif nt == "G":
                            noPAScounts[x+offset,2] += 1
                        elif nt == "T":
                            noPAScounts[x+offset,3] += 1
                        else:
                            print 'something went wrong'
                pasAssignmentFile.write("%s\t%s\n" %(cluster,"noPAS"))

    for x in range(window * 4 + 1):
        if sum(canonNTcounts[x,:]) != 0:
            canonNTcounts[x,:] = 100 * canonNTcounts[x,:] / sum(canonNTcounts[x,:]) # normalize each position to a percentage
        if sum(noncanonNTcounts[x,:]) != 0:
            noncanonNTcounts[x,:] = 100 * noncanonNTcounts[x,:] / sum(noncanonNTcounts[x,:])
        if sum(noPAScounts[x,:]) != 0:
            noPAScounts[x,:] = 100 * noPAScounts[x,:] / sum(noPAScounts[x,:])


    print canonCount
    print 100.*float(canonCount) / float(canonCount+noncanonCount+noPASCount)
    print noncanonCount
    print 100.*float(noncanonCount) / float(canonCount+noncanonCount+noPASCount)
    print noPASCount
    print 100.*float(noPASCount) / float(canonCount+noncanonCount+noPASCount)
    plt.figure(num=None,figsize=(2.777,2.5))
    plt.plot(range(-window,window+1),canonNTcounts[window:3*window+1,0],'b-',label="A")
    plt.plot(range(-window,window+1),canonNTcounts[window:3*window+1,3],'r-',label="T")
    plt.plot(range(-window,window+1),canonNTcounts[window:3*window+1,2],'g-',label="G")
    plt.plot(range(-window,window+1),canonNTcounts[window:3*window+1,1],'m-',label="C")
    plt.ylim((0,100))
    plt.ylabel("Percent bases observed")
    plt.xlabel("Relative nt position (anchored at PAS)")
    plt.title("AAUAAA")
    plt.legend()
    # plt.show()
    plt.tight_layout()
    plt.savefig(outfiles_prefix + "CanonPASprofile.pdf")
    plt.clf()
    plt.figure(num=None,figsize=(2.777,2.5))
    plt.plot(range(-window,window+1),noncanonNTcounts[window:3*window+1,0],'b-',label="A")
    plt.plot(range(-window,window+1),noncanonNTcounts[window:3*window+1,3],'r-',label="T")
    plt.plot(range(-window,window+1),noncanonNTcounts[window:3*window+1,2],'g-',label="G")
    plt.plot(range(-window,window+1),noncanonNTcounts[window:3*window+1,1],'m-',label="C")
    plt.ylim((0,100))
    plt.ylabel("Percent bases observed")
    plt.xlabel("Relative nt position (anchored at PAS)")
    plt.title("Alt PAS")
    #plt.legend()
    plt.tight_layout()
    plt.savefig(outfiles_prefix + "AltPASprofile.pdf")
    plt.clf()
    # plt.show()
    plt.clf()
    plt.figure(num=None,figsize=(2.777,2.5))
    plt.plot(range(-window,window),noPAScounts[window:3*window,0],'b-',label="A")
    plt.plot(range(-window,window),noPAScounts[window:3*window,3],'r-',label="T")
    plt.plot(range(-window,window),noPAScounts[window:3*window,2],'g-',label="G")
    plt.plot(range(-window,window),noPAScounts[window:3*window,1],'m-',label="C")
    plt.ylim((0,60))
    plt.ylabel("Percent bases observed")
    plt.xlabel("Relative nt position")
    plt.title("no PAS")
    #plt.legend()
    plt.tight_layout()
    plt.savefig(outfiles_prefix + "NoPASprofile.pdf")
    plt.clf()
    
    # plt.show()
    plt.figure(num=None,figsize=(2.777,2.5))
    plt.hist(canon_offsets,bins = [x - 0.5 for x in range(-40,max(canon_offsets))],density =True)
    plt.xlabel("Offset of PAS from called 3'UTR endpoint")
    plt.ylabel("Frequency")
    plt.title("AAUAAA")
    plt.tight_layout()
    plt.savefig(outfiles_prefix + "CanonPASoffsetsHistogram.pdf")
    # pyplot.show()
    plt.figure(num=None,figsize=(2.777,2.5))
    plt.hist(noncanon_offsets,bins = [x - 0.5 for x in range(-40,max(canon_offsets))],density = True)
    plt.xlabel("Offset of PAS from called 3'UTR endpoint")
    plt.ylabel("Frequency")
    plt.title("Alt PAS")
    plt.tight_layout()
    plt.savefig(outfiles_prefix + "AltPASoffsetsHistogram.pdf")
    # pyplot.show()


fastaFile = open(sys.argv[1],'r')
noncanonPASTable = open(sys.argv[2],'r')
pasAssignmentFile = open(sys.argv[3],'w')
outfiles_prefix = sys.argv[4]

plotAllPAS(fastaFile,noncanonPASTable,pasAssignmentFile,outfiles_prefix)
