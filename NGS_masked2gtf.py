import sys, argparse, difflib, itertools
import numpy as np
from Bio import SeqIO

# create pareser object
parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('--f', '--fasta', type=str, help='input fasta file')
parser.add_argument('--m', '--masked', type=str, help='masked fasta file')
parser.add_argument('--g', '--gtf', type=str, default="masked.gtf", help='gtf file')

# read the command line inputs
args = parser.parse_args()

# The actuall python script
#==========================================
in1 = SeqIO.parse(args.f, 'fasta')
in2 = SeqIO.parse(args.m, 'fasta')
gtf = open(args.g, 'w')

#in1 = SeqIO.parse('testB.fasta', 'fasta')
#in2 = SeqIO.parse('masked.fasta', 'fasta')


for i1,i2 in itertools.izip(in1, in2):
    old = False
    n = 0
    length = len(i1.seq)
    for n,(ii1,ii2) in enumerate(zip(str(i1.seq), str(i2.seq))):
        #n += 1
        new = int(ii1 != ii2)
        if (new == True) & (old == False):
            start = n+1
        if (new == False) & (old == True ):
            end = n
            #print start,end,i1.id
            gtf.write('%s\tmasked\texon\t%s\t%s\t.\t.\t.\tgene_id "%s"; transcript_id "%s";\n' %
                      (i1.id, start, end, i1.id ,i1.id ) )
        if (n+1 == length) & (old == True):
            end = n
            gtf.write('%s\tmasked\texon\t%s\t%s\t.\t.\t.\tgene_id "%s"; transcript_id "%s";\n' %
                      (i1.id, start, end, i1.id ,i1.id ) )           
        old = new

   

