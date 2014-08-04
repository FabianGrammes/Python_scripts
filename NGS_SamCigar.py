import argparse
import pysam
import sys

# create pareser object
parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('--s', '--sam', type=str, help='input sam file')

# read the command line inputs
args = parser.parse_args()


# The actuall python script
#==========================================

samfile = args.s

# Check for SAM file
#if samfile.split(',')[-1] != "sam":
#    print 'ERROR: script only accepts sam files'
#    sys.exit(1) # abort

def parse_cigar(Sfile):
    infile = pysam.Samfile(Sfile, "r") 
    Sout = Sfile.split(".")[0]+".out.cig.sam" 
    outfile = pysam.Samfile( Sout , "wh", template = infile )
    for e in infile: 
        if e.cigar[0][0] == 4:
            q = e.qual
            st = e.cigar[0][1]
            e.seq = e.seq[st:]
            e.qual = q[st:]
            e.cigar = e.cigar[1:]
        if e.cigar[-1][0] == 4:
            q = e.qual
            end = e.cigar[-1][1]
            e.seq = e.seq[:-end]   
            e.qual = q[:-end]
            e.cigar = e.cigar[:-1]
        outfile.write(e)

    outfile.close()

parse_cigar(samfile)
