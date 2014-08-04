import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import sys

# create pareser object
parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('--f', '--fasta', type=str, help='input fasta file')
parser.add_argument('--i', '--intersect', type=str, help='file with fasta ids to intersect; 1 id per line')
parser.add_argument('--o', '--output', type=str, default="stdout.fastq", help='output file')

# read the command line inputs
args = parser.parse_args()

# The actuall python script
#==========================================
fasta_file = args.f
int_file = args.i
out_file = args.o

ftype = fasta_file.split('.')[-1]

# read in the file with the ids to intersect
int_ids= []
int_infile = open(int_file, 'r')
for line in int_infile:
    int_ids.append(line.replace('\n',''))

int_infile.close()    

handle = open(out_file, 'w')

if ftype == 'fasta':
    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.id in int_ids:
            handle.write(">%s\n%s\n" % (record.description, record.seq))
            
    handle.close()

if ftype == 'fastq':
    for title, seq, qual in FastqGeneralIterator(open(fasta_file)):
        if title.split('/')[0] in int_ids:
            handle.write("@%s\n%s\n+\n%s\n" % (title, seq, qual))

    handle.close()
