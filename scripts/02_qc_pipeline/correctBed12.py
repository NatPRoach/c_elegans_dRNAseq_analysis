#!/usr/bin/env python2

import sys
import sets

def resolveDuplicates(str1,str2):
    ## Takes in two strings (tx_ids), and returns the root tx_id common to both.
    if len(str1) > len(str2):
        y = len(str2)
    else:
        y = len(str1)
    common_str = ""
    for x in range(y):
        if str1[x] == str2[x]:
            common_str += str1[x]
        else:
            break
    root = common_str.strip('.')
    return root
    

def getBlocks(sizes,starts):
    int_sizes = []
    int_starts = []
    for size in sizes.split(','):
        int_sizes.append(int(size))
    for start in starts.split(','):
        int_starts.append(int(start))
    return int_sizes, int_starts



def formatBlocks(int_sizes,int_starts):
    sizes = []
    starts = []
    for size in int_sizes:
        sizes.append(str(size))
    for start in int_starts:
        starts.append(str(start))
    return ",".join(sizes), ",".join(starts)



def splitAttributes(attr):
    fields = attr.split(';')
    attr_dict = {}
    for field in fields:
        key,val = field.split('=')
        attr_dict[key] = val
    return attr_dict



def testOverlap(read,single_exon_gene_list):
    putative_gene = None
    for gene in single_exon_gene_list: #gene list must be sorted
        if read[0] >= gene[0] and read[0] <= gene[1]:
            if putative_gene is not None:
                if putative_gene != gene[2]:
                    return None
            else:
                putative_gene = gene[2]
        elif read[1] >= gene[0] and read[1] <= gene[1]:
            if putative_gene is not None:
                if putative_gene != gene[2]:
                    return None
            else:
                putative_gene = gene[2]
        elif read[0] <= gene[0] and read[1] >= gene[1]:
            if putative_gene is not None:
                if putative_gene != gene[2]:
                    return None
            else:
                putative_gene = gene[2]
        elif read[1] <= gene[0]:
            return putative_gene



def getTxStructure(start,block_sizes,block_starts):
    tx_structure = []
    if len(block_sizes) > 1:
        for x in range(len(block_sizes)-1):
            size1 = block_sizes[x]
            size2 = block_sizes[x+1]
            start1 = start + block_starts[x]
            end1 = start1 + size1
            start2 = start + block_starts[x+1]
            end2 = start2 + size2
            tx_structure.append((end1,start2))
    return tuple(tx_structure)



def getTxID(chrom,strand,start,block_sizes,block_starts,tx_dict):
    tx_ids = None
    tx_structure = getTxStructure(start,block_sizes,block_starts)
    if (chrom,strand) in tx_dict:
        if tx_structure in tx_dict[(chrom,strand)]:
            tx_ids = tx_dict[(chrom,strand)][tx_structure]
    return tx_ids



def getGeneID(chrom,strand,start,block_sizes,block_starts,donor_genes,acceptor_genes):
    last_gene_id = None
    for x in range(len(block_sizes)):
        block_start = start + block_starts[x]
        block_end = block_start + block_sizes[x]
        if x == 0:
            if strand == '+':
                if (chrom,strand,block_end) in donor_genes:
                    last_gene_id = donor_genes[(chrom,strand,block_end)]
                else:
                    return None
            elif strand == '-':
                if (chrom,strand,block_end) in acceptor_genes:
                    last_gene_id = acceptor_genes[(chrom,strand,block_end)]
                else:
                    return None # Ambiguous mapping removed it
            else:
                return None #Faulty strand
        elif x == len(block_sizes) - 1:
            if strand == '+':
                if (chrom,strand,block_start) in acceptor_genes:
                    if last_gene_id != acceptor_genes[(chrom,strand,block_start)]:
                        return None # Mapping to more than one gene
                else:
                    return None # Ambiguous mapping removed it
            elif strand == '-':
                if (chrom,strand,block_start) in donor_genes:
                    if last_gene_id != donor_genes[(chrom,strand,block_start)]:
                        return None # Mapping to more than one gene
                else:
                    return None # Ambiguous mapping removed it
            else:
                return None #Faulty strand
        else:
            if strand == '+':
                if (chrom,strand,block_start) in acceptor_genes:
                    if last_gene_id != acceptor_genes[(chrom,strand,block_start)]:
                        return None # Mapping to more than one gene
                else:
                    return None # Ambiguous mapping removed it
                if (chrom,strand,block_end) in donor_genes:
                    if last_gene_id != donor_genes[(chrom,strand,block_end)]:
                        return None
                else:
                    return None
            
            elif strand == '-':
                if (chrom,strand,block_start) in donor_genes:
                    if last_gene_id != donor_genes[(chrom,strand,block_start)]:
                        return None # Mapping to more than one gene
                else:
                    return None # Ambiguous mapping removed it
                if (chrom,strand,block_end) in acceptor_genes:
                    if last_gene_id != acceptor_genes[(chrom,strand,block_end)]:
                        return None
                else:
                    return None
            else:
                return None #Faulty strand
    return last_gene_id

reference_in = open("/Users/nproach/Documents/LabFiles/Bioinformatics/NPR_Notebook/00_Data/references/ce11/c_elegans.PRJNA13758.WS265.WormBase.gff3")
window = 15

##Construct dictionary of donors and acceptor pairs
donors = {} # chr,strand -> set of donors on that chr strand
acceptors = {} #chr,strand -> set of acceptors on that chr strand
tx_id_to_gene = {}
donor_genes = {} #chr,strand,pos -> gene_id
acceptor_genes = {} #chr,strand,pos -> gene_id
tx_dict = {}
ambiguous_donors = sets.Set() #chr,strand,pos for donors with more than one gene assignment
ambiguous_acceptors = sets.Set() #chr,strand,pos for acceptors with more than one gene assignment
single_exon_genes = {}
tx_id_to_exons = {}
for line in reference_in:
    fields = line.strip().split()
    chrom = fields[0]
    if chrom == "MtDNA":
        chrom = "chrM"
    else:
        chrom = "chr" + chrom
    start = int(fields[3]) - 1
    end =  int(fields[4])
    strand = fields[6]
    if fields[2] == "exon":
        tx_id = fields[8].split(':')[1]
        index = tx_id.find('.') + 1
        if tx_id in tx_id_to_exons:
            tx_id_to_exons[tx_id].append((chrom,strand,start,end))
        else:
            tx_id_to_exons[tx_id] = [(chrom,strand,start,end)]
    elif fields[2] == "gene":
        attr = splitAttributes(fields[8])
        if "Name" in attr:
            gene_id = attr["Name"]
            if "sequence_name" in attr:
                tx_id = attr["sequence_name"]
                tx_id_to_gene[tx_id] = gene_id
            else:
                print "something went wrong"
        else:
            print "something went wrong"



for tx_id in tx_id_to_exons: ##WARNING: GFF3 file must be sorted for this approach to work
    exons = tx_id_to_exons[tx_id]
    index = tx_id.find('.') + 1
    t_flag = True
    while True:
        try:
            if t_flag: #For some reason some genes have a t instead of a # immediately after the 1st '.' because this format sucks
                if tx_id[index] != 't': #Need to handle that, so this is my solution.
                    int(tx_id[index])
            else:
                int(tx_id[index])
            index += 1
            t_flag = False
        except:
            break
    tx_id2 = tx_id[:index]
    gene_id = tx_id_to_gene[tx_id2]
    if len(exons) > 1:
        chrom = exons[0][0]
        strand = exons[0][1]
        intron_structure = []
        for x in range(len(exons)-1): ##ignore 5' truncations for now.
            start = exons[x][3]
            end = exons[x+1][2]
            intron_structure.append((start,end))
        if (chrom,strand) in tx_dict:
            if tuple(intron_structure) in tx_dict[(chrom,strand)]:
                # reduced_tx_id = resolveDuplicates(tx_id,tx_dict[(chrom,strand)][tuple(intron_structure)][0])
                tx_dict[(chrom,strand)][tuple(intron_structure)].append(tx_id)
                # tx_dict[(chrom,strand)][tuple(intron_structure)] = [reduced_tx_id]
            else:
                tx_dict[(chrom,strand)][tuple(intron_structure)] = [tx_id]
        else:
            tx_dict[(chrom,strand)]= {tuple(intron_structure) : [tx_id]}
            
        # for x in range(len(intron_structure)):
        #     if strand == '+':
        #         temp_intron_structure = intron_structure[x:]
        #     elif strand == '-':
        #         if x == 0:
        #             temp_intron_structure = intron_structure[:]
        #         else:
        #             temp_intron_structure = intron_structure[:-x]
        #     if (chrom,strand) in tx_dict:
        #         if tuple(temp_intron_structure) in tx_dict[(chrom,strand)]:
        #             tx_dict[(chrom,strand)][tuple(temp_intron_structure)].append(tx_id)
        #         else:
        #             tx_dict[(chrom,strand)][tuple(temp_intron_structure)] = [tx_id]
        #     else:
        #         tx_dict[(chrom,strand)]= {tuple(temp_intron_structure) : [tx_id]}
        for i, exon in enumerate(exons):
            chrom,strand,start,end = exon
            if (chrom,strand) not in donors:
                donors[(chrom,strand)] = sets.Set()
                acceptors[(chrom,strand)] = sets.Set()
            if i == 0: #Only add the donor for positive strand, acceptor for negative
                if strand == '+':
                    donors[(chrom,strand)].add(end)
                    if (chrom,strand,end) not in donor_genes:
                        donor_genes[(chrom,strand,end)] = gene_id#chr,strand,pos -> gene_id
                    elif gene_id != donor_genes[(chrom,strand,end)]:
                        ambiguous_donors.add((chrom,strand,end))
                elif strand == '-':
                    acceptors[(chrom,strand)].add(end)
                    if (chrom,strand,end) not in acceptor_genes:
                        acceptor_genes[(chrom,strand,end)] = gene_id #chr,strand,pos -> gene_id
                    elif gene_id != acceptor_genes[(chrom,strand,end)]:
                        ambiguous_acceptors.add((chrom,strand,end))
                else:
                    print "Something went wrong, strand not + or -"
            elif i == len(exons) - 1: # Only add the acceptor for + strand, donor for negative
                if strand == '+':
                    acceptors[(chrom,strand)].add(start)
                    if (chrom,strand,start) not in acceptor_genes:    
                        acceptor_genes[(chrom,strand,start)] = gene_id #chr,strand,pos -> gene_id
                    elif gene_id != acceptor_genes[(chrom,strand,start)]:
                        ambiguous_acceptors.add((chrom,strand,start))
                elif strand == '-':
                    donors[(chrom,strand)].add(start)
                    if (chrom,strand,start) not in donor_genes:
                        donor_genes[(chrom,strand,start)] = gene_id#chr,strand,pos -> gene_id
                    elif gene_id != donor_genes[(chrom,strand,start)]:
                        ambiguous_donors.add((chrom,strand,start))
                else:
                    print "Something went wrong, strand not + or -"
            else: #Add both donor and acceptor
                if strand == '+':
                    acceptors[(chrom,strand)].add(start)
                    if (chrom,strand,start) not in acceptor_genes:    
                        acceptor_genes[(chrom,strand,start)] = gene_id #chr,strand,pos -> gene_id
                    elif gene_id != acceptor_genes[(chrom,strand,start)]:
                        ambiguous_acceptors.add((chrom,strand,start))
                    donors[(chrom,strand)].add(end)
                    if (chrom,strand,end) not in donor_genes:
                        donor_genes[(chrom,strand,end)] = gene_id#chr,strand,pos -> gene_id
                    elif gene_id != donor_genes[(chrom,strand,end)]:
                        ambiguous_donors.add((chrom,strand,end))
                elif strand == '-':
                    acceptors[(chrom,strand)].add(end)
                    if (chrom,strand,end) not in acceptor_genes:
                        acceptor_genes[(chrom,strand,end)] = gene_id #chr,strand,pos -> gene_id
                    elif gene_id != acceptor_genes[(chrom,strand,end)]:
                        ambiguous_acceptors.add((chrom,strand,end))
                    donors[(chrom,strand)].add(start)
                    if (chrom,strand,start) not in donor_genes:
                        donor_genes[(chrom,strand,start)] = gene_id#chr,strand,pos -> gene_id
                    elif gene_id != donor_genes[(chrom,strand,start)]:
                        ambiguous_donors.add((chrom,strand,start))
                else:
                    print "Something went wrong"
    else: #Ideally we'd want a way to match single exon reads with single exon genes
        chrom,strand,start,end = exons[0]
        if (chrom,strand) in single_exon_genes:
            single_exon_genes[(chrom,strand)].append((start,end,gene_id))
        else:
            single_exon_genes[(chrom,strand)] = [(start,end,gene_id)]
    
    for key in single_exon_genes:
        single_exon_genes[key].sort()


for donor in ambiguous_donors:
    del donor_genes[donor]
for acceptor in ambiguous_acceptors:
    del acceptor_genes[acceptor]


in_file = open(sys.argv[1])
edit_file = open(sys.argv[2] + ".edit.txt",'w')
removed_file = open(sys.argv[2] + ".removed.txt",'w')
gene_file = open(sys.argv[2] + ".gene.txt",'w')

tmp1 = [-x for x in range(1,window+1)]
tmp2 = [ x for x in range(1,window+1)]
window_range = []
for x in range(len(tmp1)):
    window_range.append(tmp1[x])
    window_range.append(tmp2[x])
for line in in_file:
    fields= line.strip().split()
    if len(fields) == 6: ##If bed6 aka no introns, it passes by definition
        print line.strip()
        continue
    elif fields[9] == '1': ## Again if no introns, passes by definition
        chrom = fields[0]
        start = int(fields[1])
        strand = fields[5]
        block_sizes,block_starts = getBlocks(fields[10],fields[11]) 
        gene_id = testOverlap((start,start+block_sizes[0]), single_exon_genes[(chrom,strand)])
        if gene_id is not None:
            gene_file.write("%s\t%s\n" %(fields[3],gene_id))
            print line.strip()
        continue
        
    chrom = fields[0]
    start = int(fields[1])
    end =  int(fields[2])
    block_sizes,block_starts = getBlocks(fields[10],fields[11]) 
    strand = fields[5]
    keep_flag = True
    edit_flag = False
    donor_edits = []
    acceptor_edits = []
    if (chrom,strand) not in acceptors: ## shouldn't ever be the case, but worth testing for
        #print "Something went wrong"
        #print chrom
        #print strand
        continue
    if strand == "+":
        for x in range(1,len(block_starts)):
            if start + block_starts[x] not in acceptors[(chrom,strand)]:
                keep_flag = False
                correction = None
                #for y in (range(-window,0) + range(1,window+1)):
                for y in window_range:
                    z = y + start + block_starts[x]
                    if z in acceptors[(chrom,strand)]:
                        keep_flag = True
                        correction = y
                        break
                if keep_flag:
                    if block_starts[x] + correction >= 0 and block_sizes[x] - correction > 0: ## Possible source of error
                        block_starts[x] = block_starts[x] + correction
                        block_sizes[x] = block_sizes[x] - correction
                        if correction != 0:
                            edit_flag = True
                            acceptor_edits.append(correction)
                    else:
                        keep_flag = False
                        break
                else:
                    break
            if not keep_flag:
                break
            if start+block_starts[x-1] + block_sizes[x-1] not in donors[(chrom,strand)]:
                keep_flag = False
                correction = None
                #for y in (range(-window,0) + range(1,window+1)):
                for y in window_range:
                    z = y + start + block_starts[x-1] + block_sizes[x-1]
                    if z in donors[(chrom,strand)]:
                        keep_flag = True
                        correction = y
                        break
                if keep_flag:
                    if block_sizes[x-1] + correction > 0:
                        block_sizes[x-1] = block_sizes[x-1] + correction
                        if correction != 0:
                            edit_flag = True
                            donor_edits.append(correction)
                    else:
                        keep_flag = False
                        break ## Possible source of error??
                else:
                    break
    elif strand == "-":
        for x in range(1,len(block_starts)):
            if start + block_starts[x] not in donors[(chrom,strand)]:
                keep_flag = False
                correction = None
                #for y in (range(-window,0) + range(1,window+1)):
                for y in window_range:
                    z = y + start+block_starts[x]
                    if z in donors[(chrom,strand)]:
                        keep_flag = True
                        correction = y
                        break
                # if fields[3] == "068acb24-3e36-4106-9eb8-5ad315dab20a":
                #     print keep_flag
                if keep_flag:
                    if block_starts[x] + correction >= 0 and block_sizes[x] - correction > 0: ## Possible source of error
                        block_starts[x] = block_starts[x] + correction
                        block_sizes[x] = block_sizes[x] - correction
                        if correction != 0:
                            edit_flag = True
                            donor_edits.append(-1*correction)
                    else:
                        keep_flag = False
                        break
                else:
                    break
            if start+block_starts[x-1] + block_sizes[x-1] not in acceptors[(chrom,strand)]:
                keep_flag = False
                correction = None
                #for y in (range(-window,0) + range(1,window+1)):
                for y in window_range:
                    z = y + start + block_starts[x-1] + block_sizes[x-1]
                    if z in acceptors[(chrom,strand)]:
                        keep_flag = True
                        correction = y
                        break
                if keep_flag:
                    if block_sizes[x-1] + correction > 0:
                        block_sizes[x-1] = block_sizes[x-1] + correction
                        if correction != 0:
                            edit_flag = True
                            acceptor_edits.append(-1*correction)
                    else:
                        keep_flag = False
                        break
                else:
                    break
    gene_id = getGeneID(chrom,strand,start,block_sizes,block_starts,donor_genes,acceptor_genes)
    tx_ids = getTxID(chrom,strand,start,block_sizes,block_starts,tx_dict)
    if gene_id is not None:
        # print gene_id
        if tx_ids is None:
            tx_ids = []
        if keep_flag:
            fields[10], fields[11] = formatBlocks(block_sizes,block_starts)
            print "\t".join(fields)
            if len(tx_ids) == 0:
                common_tx_id = ""
            else:
                common_tx_id = tx_ids[0]
            if len(tx_ids) > 1:
                for x in range(1,len(tx_ids)):
                    common_tx_id = resolveDuplicates(common_tx_id,tx_ids[x])
                
            gene_file.write("%s\t%s\t%s\t%s\n" %(fields[3],gene_id,",".join([str(x) for x in list(tx_ids)]), common_tx_id))
            if edit_flag:
                if len(donor_edits) != 0:
                    if len(acceptor_edits) != 0:
                        edit_file.write("%s\t%s\t%s\n" %(fields[3],','.join(str(i) for i in donor_edits),','.join(str(i) for i in acceptor_edits)))
                    else:
                        edit_file.write("%s\t%s\t%s\n" %(fields[3],','.join(str(i) for i in donor_edits),'.'))
                else:
                    edit_file.write("%s\t%s\t%s\n" %(fields[3],'.',','.join(str(i) for i in acceptor_edits)))
        else:
            removed_file.write("%s\n" %(fields[3]))
    else:
        removed_file.write("%s\n" %(fields[3]))