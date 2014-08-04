import os

'''
Function to parse fpkm_tracking files produced
by cufflinks. The function will reaturn the number
of genes above a certain treshold -- binning...
''' 

def parse_fpkm(FILE):
    c1 = c2 = c3 = c4 = 0
    val =[]
    infile = open(FILE, 'r')
    infile.next()
    for line in infile:
        fpkm = float(line.split('\t')[9])
        if fpkm > 100:
            c4 += 1
        if fpkm > 10:
            c3 += 1 
        if fpkm > 1:
            c2 += 1
        if fpkm > 0:
            c1 += 1     
    result = [FILE,c1,c2,c3,c4]
    return result
    infile.close()

# list files
list_files = []
for track in os.listdir(os.getcwd()):
    if track.endswith(".fpkm_tracking"):
        list_files.append(track)

outfile = open('Binned_FPKM_Results.tab' , "w")
outfile.write('%s\t%s\t%s\t%s\t%s\n' % ('FILE', 'FPKM>0','FPKM>1','FPKM>10','FPKM>100'))
for ent in list_files:
    a = parse_fpkm(ent)
    outfile.write('%s\t%s\t%s\t%s\t%s\n' % (a[0], a[1], a[2], a[3], a[4]))
    
outfile.close() 





            
        

        
