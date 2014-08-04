import argparse
import pysam
import sys

# create pareser object
parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('--m', '--sam_m', type=str, help='sam with mapped sequences')
parser.add_argument('--u', '--sam_u', type=str, help='sam with unmapped sequences (No SAM header needed)')
parser.add_argument('--o', '--out', type=str, default = "itersesect.fastq", help='outfile')

# read the command line inputs
args = parser.parse_args()

'''
The script takes 2 SAM files and intersects the reads that can be
found in both of them.

The second SAM file should contain only alignment sections from the
SAM file e.g the file should be generated by something like:
$ samtools view -S -F4 sample.sam > sample.unmapped.sam

The output is in FASTA format.
'''



# The actuall python script
#=======================================
map_SAM = args.m
unmap_SAM = args.u
out_SAM = args.o

sam_unmap = []

u_file = open(unmap_SAM, 'r')
for line in u_file:
    sam_unmap.append(line.split('\t')[0])

u_file.close()

count = 1

m_file = pysam.Samfile(map_SAM, "r")
outfile = open(out_SAM , "w") 
for entry in m_file:
    if entry.qname in sam_unmap:
        count += 1
        outfile.write('>%s<==>%s\n%s\n' % (entry.qname, m_file.getrname(entry.rname), entry.seq))  
        if count > 20:
            break
        
outfile.close()

