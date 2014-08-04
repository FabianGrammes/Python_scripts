import argparse
from Bio import SeqIO
from Bio.Seq import Seq

'''
This script does 2 things:
1: It removes empty reacords from FASTA
2: It concatenates the FASTA name 
'''

parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('--f', '--fasta', type=str, help='input fasta file')
parser.add_argument('--o', '--out', type=str, help='output fasta file')

# read the command line inputs
args = parser.parse_args()

# The actuall python script
#==========================================
empty=0
rec=0
out=[]

for record in SeqIO.parse(args.f, "fasta"):
    record.id = '|'.join(record.description.split())
    record.description = ""
    if str(record.seq) == "":
        empty+=1 
    else:
        out.append(record)
        rec+=1
             
SeqIO.write(out, args.o, "fasta")

