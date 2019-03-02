#!/usr/bin/env python2

import sys

## For now assume bed12
infile = open(sys.argv[1])
## Outputting bed6 in v2
outfile = open(sys.argv[2],'w')
window = 60 

for line in infile:
    fields = line.strip().split()
    start = int(fields[1])
    end  = int(fields[2])
    strand = fields[5]
    # block_sizes = [int(x) for x in fields[10].split(',')]
    # block_starts = [int(x) for x in fields[11].split(',')]
    if end - start < window:
        difference = window - (end - start)
        if strand == '+':
            start -= difference
            end += window
        elif strand == '-':
            end += difference
            start -= window
    else:
        if strand == '+':
            end += window
        elif strand == '-':
            start -= window
    start = str(start)
    end = str(end)
    line_out = "%s\t%s\t%s\t%s\t%s\t%s\n" %(fields[0],start,end,fields[3],fields[4],fields[5])
    outfile.write(line_out)

# for line in infile:
#     fields = line.strip().split()
#     start = int(fields[1])
#     end  = int(fields[2])
#     strand = fields[5]
#     block_sizes = [int(x) for x in fields[10].split(',')]
#     block_starts = [int(x) for x in fields[11].split(',')]
#     if sum(block_sizes) < window:
#         difference = window - sum(block_sizes)
#         if strand == '+':
#             start -= difference
#             block_sizes[0] += difference
#             new_block_starts = [0]
#             for x in block_starts[1:]:
#                 new_block_starts.append(x+difference)
#             block_starts = new_block_starts
#
#             end += window
#             block_sizes[-1] += window
#
#         elif strand == '-':
#             end += difference
#             block_sizes[-1] += difference
#
#             start -= window
#             block_sizes[0] += window
#             new_block_starts = [0]
#             for x in block_starts[1:]:
#                 new_block_starts.append(x+window)
#             block_starts = new_block_starts
#     else:
#         if strand == '+':
#             end += window
#             block_sizes[-1] += window
#
#         elif strand == '-':
#             start -= window
#             block_sizes[0] += window
#             new_block_starts = [0]
#             for x in block_starts[1:]:
#                 new_block_starts.append(x+window)
#             block_starts = new_block_starts
#
#     start = str(start)
#     end = str(end)
#     block_sizes = ','.join([str(x) for x in block_sizes])
#     block_starts = ','.join([str(x) for x in block_starts])
#     line_out = "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %(fields[0],start,end,fields[3],fields[4],fields[5],start,start,fields[8],fields[9],block_sizes,block_starts)
#     outfile.write(line_out)