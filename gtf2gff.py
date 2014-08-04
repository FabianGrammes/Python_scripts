import HTSeq,argparse

"""
The following script will reformat a GTF file (assumed to come from cufflinks)
to GFF format.

The script also makes it possible to intersect entries from a ceratin locus e.g
gtf:gene_id or gtf:transcript_id
"""
#=====================================================================

# create pareser object
parser = argparse.ArgumentParser()

# add command line options
parser.add_argument('--i', '--gtf', type=str, help='input .GTF file')
parser.add_argument('--o', '--out', type=str, default='out.gff', help='output .GFF file')
parser.add_argument('--r', '--gene', type=str, default='empty', help='Gene to intersect')

# read the command line inputs
args = parser.parse_args()


#=====================================================================
class gtf2gff:
    def __init__(self, infile, outfile, gene):
        self.infile = infile
        self.out = open(outfile, 'w')
        if gene == 'empty':          
            self.__parser()
        elif gene != 'empty':
            self.__parserG(gene)
        self.out.close()

    def __write_extra(self, line, start, end, etype):
        '''
        function to write GFF line for genes or transcripts(mRNA)
        '''
        if etype == "gene":
            attr = 'ID='+line.name+';Name='+line.name
        if etype == "mRNA":
            attr = 'ID='+line.attr['transcript_id']+';Parent='+line.name
        self.out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(line.iv.chrom, "py-gtf2gff", etype, start, end,".", line.iv.strand,".", attr))

    def __write_entry(self, line):
        '''
        function to write GFF line for exons
        '''
        attr = 'ID='+line.attr['transcript_id']+'.'+line.attr['exon_number']+';Parent='+line.attr['transcript_id']
        self.out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' %(line.iv.chrom, "py-gtf2gff", line.type, line.iv.start+1, line.iv.end, "." , line.iv.strand,".", attr))

    def __transcript_handle(self, gene_box):
        '''
        subroutine collecting all exons for a transcript and writing them
        to file 
        '''
        mRNA_start = []
        mRNA_end = []
        mRNA_box = []
        old_entry = ''
        for entry in gene_box:
            if old_entry == '':
                pass
            elif entry.attr['transcript_id'] == old_entry.attr['transcript_id']:
                mRNA_start.append(old_entry.iv.start+1)
                mRNA_end.append(old_entry.iv.end)
                mRNA_box.append(old_entry)
            elif entry.attr['transcript_id'] != old_entry.attr['transcript_id']:
                mRNA_start.append(old_entry.iv.start+1)
                mRNA_end.append(old_entry.iv.end)
                mRNA_box.append(old_entry)
                self.__write_extra( old_entry, min(mRNA_start), max(mRNA_end), 'mRNA' )
                for exon in mRNA_box:
                   self.__write_entry( exon)
                mRNA_start = []
                mRNA_end = []
                mRNA_box = []
            old_entry = entry
        mRNA_start.append(entry.iv.start+1)
        mRNA_end.append(entry.iv.end)
        mRNA_box.append(entry)
        self.__write_extra( entry, min(mRNA_start), max(mRNA_end), 'mRNA' )
        for exon in mRNA_box:
            self.__write_entry(exon )  


    def __parser(self):
        '''
        FINAL FUNCTION!
        parses the gtf and uses the functions from the class to make
        a GFF file
        '''
        gtf_file = HTSeq.GFF_Reader(self.infile, end_included=True )
        gene_start = []
        gene_end = []
        gene_box = []
        old_line = ''
        for line in gtf_file:
            if old_line == '':
                pass
            elif old_line.name == line.name:
                gene_box.append(old_line)
                gene_start.append(old_line.iv.start+1)
                gene_end.append(old_line.iv.end)
            elif old_line.name != line.name:
                gene_box.append(old_line)
                gene_start.append(old_line.iv.start+1)
                gene_end.append(old_line.iv.end)
                self.__write_extra( old_line, min(gene_start), max(gene_end), 'gene')
                self.__transcript_handle( gene_box)
                gene_start = []
                gene_end = []
                gene_box = []
            old_line = line
        gene_box.append(line)
        gene_start.append(line.iv.start+1)
        gene_end.append(line.iv.end)
        self.__write_extra( line, min(gene_start), max(gene_end), 'gene')
        self.__transcript_handle( gene_box)

    def __parserG(self, gene):
        '''
        FINAL FUNCTION!Intersect gene
        parses the gtf and uses the functions from the class to make
        a GFF file
        '''
        gtf_file = HTSeq.GFF_Reader(self.infile, end_included=True )
        gene_start = []
        gene_end = []
        gene_box = []
        old_line = ''
        for line in gtf_file:
            if old_line == '':
                pass
            elif (old_line.name == line.name) :
                gene_box.append(old_line)
                gene_start.append(old_line.iv.start+1)
                gene_end.append(old_line.iv.end)
            elif (old_line.name != line.name):
                gene_box.append(old_line)
                gene_start.append(old_line.iv.start+1)
                gene_end.append(old_line.iv.end)
                if(old_line.name == gene):
                    self.__write_extra( old_line, min(gene_start), max(gene_end), 'gene')
                    self.__transcript_handle( gene_box)
                gene_start = []
                gene_end = []
                gene_box = []
            old_line = line
        if (line.name == gene):
            gene_box.append(line)
            gene_start.append(line.iv.start+1)
            gene_end.append(line.iv.end)
            self.__write_extra( line, min(gene_start), max(gene_end), 'gene')
            self.__transcript_handle( gene_box)
            
"===================================================================="

# Call the class

gtf2gff(args.i, args.o, args.r)

