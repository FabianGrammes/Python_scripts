import pybedtools,os,HTSeq

'''

'''
"============= Functions ============================================="

class gtf:
    def __init__(self, file_list):
        '''
        Function counts the intersects between for all gtf files supplied
        by the file_list object. Intersect is only accepted when the strand
        information matches. 
        '''
        self.file_list = file_list
        row_val = []
        for r in self.file_list:
            col_val = []
            r_file = pybedtools.BedTool(r)
            for c in file_list:
                c_file = pybedtools.BedTool(c)
                if r_file == c_file:
                    nr = r_file.count()
                else:
                    nr = r_file.intersect(c_file, s=True, v=True).count()
                col_val.append(nr)
            row_val.append(col_val)
        self.matrix = row_val
                
    def printout(self, outfile):
        '''
        Fuction to print the intersection matrix to file
        '''
        outfile = open(outfile, 'w' )
        outfile.write('Row-in-Col\t%s\n' % ('\t'.join(map(str, self.file_list))))
        for i in range(len(self.file_list)):
            outfile.write('%s\t%s\n' % (self.file_list[i], '\t'.join(map(str, self.matrix[i]))))
        outfile.close()


class gtfparse:
    def __init__(self, infile):
        self.infile = infile
        self.__parser()
        
    def __write(self, gff_line, out, start , end):
        '''
        write gtf file function...
        '''
        out.write('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\tgene_id "%s"; Nr.exons %s\n' % (gff_line.iv.chrom, "py-gtf", "gene", min(start), max(end),".", gff_line.iv.strand, ".", gff_line.name, len(start)))

    def __parser(self):
        '''
        Function generates a .gtf file containing the start and end position
        of a gene, defined by gene_id; the function also returns the total
        number of exons per gene.
        The file gets the extension *.genes.gtf. Normmally .gtf generated
        via cuffmerge are used as input. 
        '''
        out = self.infile.split(".")[0]+".genes.gtf"
        outfile = open(out, 'w')
        gtf_file = HTSeq.GFF_Reader(self.infile, end_included=True )
        st_list = []
        end_list = []
        trans =[]
        old =''
        line_nr=0
        for line in gtf_file:
            line_nr+=1
            if old == '':
                pass
            elif old.name == line.name:
                st_list.append(old.iv.start+1)
                end_list.append(old.iv.end)
            elif old.name != line.name:
                st_list.append(old.iv.start+1)
                end_list.append(old.iv.end)
                try:
                    self.__write(old, outfile, st_list, end_list)
                    st_list = []
                    end_list = []
                except:
                    print 'ERROR AT:\t%s\tLine Nr:%d' % (old.iv.chrom, line_nr)
            old = line
        if old.name != line.name:
            self.__write(line, outfile, line.iv.start+1, line.iv.end )
        elif old.name == line.name:
            st_list.append(line.iv.start+1)
            end_list.append(line.iv.end)
            self.__write(line, outfile, st_list, end_list )
            outfile.close()

    
"===================================================================="

'''
Steps to analyze the gtfs:
1: get all .gtf files
2: get the intersection matrix for all exons (using stranded info)
3: Use all exons from all transcripts to get the genomic location
   for each gene
4: get intersection for gene genomic locations (using stranded info)
'''
file_list = []
for i in os.listdir( os.getcwd()):
    if i.split('.')[-1] == 'gtf':
        file_list.append(i)

exons = gtf(file_list)
exons.printout('exon_intersect_matrix.txt')

for entry in file_list:
    gtfparse(entry)

file_list2 = []
for i in os.listdir( os.getcwd()):
    if i.split('.')[-2] == 'genes':
        file_list2.append(i)

genes = gtf(file_list2)
genes.printout('genes_intersect_matrix.txt')

