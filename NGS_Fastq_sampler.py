import sys, random, itertools, argparse, gzip
from Bio import SeqIO


# create pareser object
parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('--f1', '--fastq_1', type=str, help='Read 1')
parser.add_argument('--f2', '--fastq_2', type=str, help='Read 2')
parser.add_argument('--p', '--prop', type=float, default=0.1, help='proportion to sample')
parser.add_argument('--o1', '--out1', type=str, default="OUT_1.fastq", help='outputfile 1')
parser.add_argument('--o2', '--out2', type=str, default="OUT_2.fastq", help='outputfile 2')

# read the command line inputs
args = parser.parse_args()

'''
The script takes 2 FASTQ files (read1 and 2- same order);
and takes a random subset of size p and writes them
in FASTA format to 1 File
'''

# The actuall python script
#==========================================


in1 = SeqIO.parse( gzip.open(args.f1, 'r'), 'fastq')
in2 = SeqIO.parse( gzip.open(args.f2, 'r'), 'fastq')
out1 = open(args.o1, 'w')
out2 = open(args.o2, 'w')

fraction = args.p

for read1, read2 in itertools.izip( in1, in2 ):
    if random.random() < fraction:
        out1.write( read1.format("fastq") )
        out2.write( read2.format("fastq") )

out1.close()
out2.close()
