import argparse
from Bio import SearchIO

"""
Script to parse a BLAST+ .xml file into a tab file.
Advatage is that you can get the description information
on query/hit, which is not possible otherwise.

OUTPUT with header
"""

#=====================================================================

# create pareser object
parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('--i', '--xml', type=str, help='input blast .xml file')
parser.add_argument('--o', '--out', type=str, default='blast_hits.tab', help='output .tab file')

# read the command line inputs
args = parser.parse_args()

#=====================================================================
outfile = open(args.o, "w")

outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                          ('q.id', 'q.description', 's.id', 's.description', 'q.start',
                           'q.end', 'q.length', 's.start', 's.end', 's.length',
                           'evalue', 'align_length'))

for qid in SearchIO.parse(args.i, 'blast-xml'):
    for hit in qid:
        for hsp in hit:
            outfile.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" %
                          (qid.id, qid.description, hit.id, hit.description, hsp.query_start,
                           hsp.query_end, qid.seq_len, hsp.hit_start, hsp.hit_end, hit.seq_len,
                           hsp.evalue, hsp.aln_span))

outfile.close()













