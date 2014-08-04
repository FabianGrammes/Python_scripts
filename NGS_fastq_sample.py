import sys, random, itertools, argparse
from Bio import SeqIO


# create pareser object
parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('--f1', '--fastq_1', type=str, help='Read 1')
parser.add_argument('--f2', '--fastq_2', type=str, help='Read 2')
parser.add_argument('--p', '--prop', type=float, default=0.1, help='proportion to sample')
parser.add_argument('--o', '--out', type=str, default="OUT.fatsa", help='outputfile')

# read the command line inputs
args = parser.parse_args()

'''
The script takes 2 FASTQ files (read1 and 2- same order);
and takes a random subset of size p and writes them
in FASTA format to 1 File
'''

# The actuall python script
#==========================================

in1 = SeqIO.parse(args.f1, 'fastq')
in2 = SeqIO.parse(args.f2, 'fastq')
out = open(args.o, 'w')

fraction = args.p

for read1, read2 in itertools.izip( in1, in2 ):
    if random.random() < fraction:
        out.write(">%s\n%s\n" % (read1.id.split('\t')[0], read1.seq))
        out.write(">%s\n%s\n" % (read2.id.split('\t')[0], read2.seq))

out.close()
