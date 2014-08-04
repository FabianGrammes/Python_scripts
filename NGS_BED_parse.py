import argparse
from Bio import SeqIO
from Bio.Seq import Seq
import sys

# create pareser object
parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('--b', '--fasta', type=str, help='input .bed file')
parser.add_argument('--s', '--strand', type=str, help='if only entries from the +/- strand should be retained')
parser.add_argument('--o', '--output', type=str, help='output file')


# read the command line inputs
args = parser.parse_args()

# The actuall python script
#==========================================
bed_file = args.b
strand = args.s
out_file = args.o




handle = open(out_file, 'w')
handle.write("%s\t%s\t%s\t%s\t%s\t%s\n" % ("transcript_id", "transcript_length", "orf_start", "orf_end", "length_orf", "strand"))

int_infile = open(bed_file, 'r')
int_infile.next() # skip the 1st line
for line in int_infile:
    record = line.split('\t')
    aa_hit = record[3].split('|')[-1].split(':')[-1].split('_')[0]
    if record[5] == strand:
        handle.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (record[0], record[2], record[6], record[7], aa_hit, record[5]))

int_infile.close()    
handle.close()    

