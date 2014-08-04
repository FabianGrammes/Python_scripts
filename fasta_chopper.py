import os
import argparse
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import sys

# create pareser object
parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('--f', '--fasta', type=str, help='input fasta file')
parser.add_argument('--n', '--not', type=str, help='FASTA Ids of the files to be removed')

# read the command line inputs
args = parser.parse_args()

# The actuall python script
#==========================================
fasta_file = args.f
remove_ids = args.n
#chunks = args.c


def write_fasta(dir,record, Nr):
    file_path = dir+record.name+".fasta"
    fileconn = open(file_path, 'w')
    fileconn.write(">%s\n%s\n" % (record.description, record.seq))
    fileconn.close()    

os.makedirs("Chunk1")
os.makedirs("Chunk2")
os.makedirs("Chunk3")
os.makedirs("Chunk4")

# read in the file with the ids that should be removed
rem_ids= []
rem_infile = open(remove_ids, 'r')
for line in rem_infile:
  rem_ids.append(line.replace('\n',''))

r=0 
c=0
for record in SeqIO.parse(fasta_file, "fasta"):
    c+=1
    if record.name in rem_ids:
        r+=1
        pass
    else:
        if c <= 16000:
            write_fasta("Chunk1/", record, c)
        elif c <= 32000:
            write_fasta("Chunk2/", record, c)
        elif c <= 48000:
            write_fasta("Chunk3/", record, c)
        elif c <= 70000:
            write_fasta("Chunk4/", record, c)

print "processed %d FASTA entries; %d entries removed" % (c,r)
