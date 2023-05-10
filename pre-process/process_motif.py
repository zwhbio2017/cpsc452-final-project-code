# Usage: Select motifs that its TF is in the gene list of our input data and output them to another file

import sys

input_motif = sys.argv[1]
input_TSS = sys.argv[2]
output_motif = sys.argv[3]
input_bed = sys.argv[4]
output_bed = sys.argv[5]

f1 = open(input_motif, 'r')
lines_motif = f1.readlines()
list_motif = [line.strip() for line in lines_motif]

f2 = open(input_TSS, 'r')
lines_TSS = f2.readlines()
list_gene = list(set([line.strip().split('\t')[0] for line in lines_TSS[1:]]))

d_motif = {}

for motif in list_motif:
    if motif.split('(')[0].split(':')[0].upper() in list_gene and motif not in d_motif.keys():
        d_motif[motif] = motif.split('(')[0].split(':')[0].upper()

output_content = []

for motif in d_motif.keys():
    output_content.append(motif + '\t' + d_motif[motif] + '\n')
    
f3 = open(output_motif, 'w')
f3.writelines(output_content)
f3.close()

output_content = []

with open(input_bed, 'r') as f4:
    line = f4.readline()
    i = 1
    while line:
        line = line.strip().split('\t')
        if line[3] in d_motif.keys():
            line.append(d_motif[line[3]])
            output_content.append('\t'.join(line) + '\n')
        if len(output_content) >= 1000000:
            print('Processed ' + str(i) + ' lines.')
            f5 = open(output_bed, 'a')
            f5.writelines(output_content)
            f5.close()
            output_content = []
        i += 1
        line = f4.readline()

