import sys
from collections import defaultdict

'''Get counts per contig from sam file'''

input_sam = sys.argv[1]
output_counts = input_sam.rsplit('.', 1)[0] + '.counts.tsv'

d = defaultdict(int)
with open(input_sam) as sam, open(output_counts, 'w') as counts:
    for line in sam:
        if line[0] == '@':
            continue
        contig = line.split('\t')[2]
        d[contig] += 1
    
    for contig in d:
        counts.write(f'{contig}\t{d[contig]}\n')
