import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import sys

# create pareser object
parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('--f', '--fasta', type=str, help='input fasta file')
parser.add_argument('--o', '--output', type=str, default='empty', help='output file')

# read the command line inputs
args = parser.parse_args()

# The actuall python script
#==========================================
fasta_file = args.f
out_file = args.o

if out_file == 'empty':
    lv = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        lv.append(len(record.seq))

    print """
             Nr. of sequences: %d
             Min:              %d
             1st Quartil:      %d 
             Mean length:      %.1f 
             Median length:    %.1f
             3rd Quartil:      %d 
             Max:              %d 
             """ % (len(lv), np.min(lv),np.percentile(lv, 25),  np.mean(lv), np.median(lv), np.percentile(lv,75), np.max(lv))

if out_file != 'empty':
    handle = open(out_file, 'w')
    handle.write("%s\t%s\t%s\n" % ("transcript_id", "gene_id", "transcript_length"))

    for record in SeqIO.parse(fasta_file, "fasta"):
        id_2 = record.description.split(' ')[1].replace('gene=', '')
        length = len(record.seq)
        handle.write("%s\t%s\t%d\n" % (record.id, id_2, length))

## handle.close()    
                                        
    
     
