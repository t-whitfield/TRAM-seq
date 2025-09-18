from optparse import OptionParser

# Report only for appropriate bases.  Skip others.
# Note that the number of skips will be strand-specific, leading to
# pileup files of different lengths on the positive and negative strands.

# Translate bcftools output into custom pileup format with per-nucleotide
# mismatch and deletion ratios reported.

# Troy Whitfield, 27 Jan. 2022.

def basetab(x,y):

    bases={'A':0,'C':0,'G':0,'T':0}
    z = list(bases.keys())
    for i in range(0,len(z)):
        if z[i] in x:
            bases[z[i]]=y[x.index(z[i])]
        
    return bases

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("--pileup", dest="PILEUP", help="input taken from bfctools mpileup custom output")
    parser.add_option("--outfile1", dest="OUTFILE1", help="output file/path")
    parser.add_option("--outfile2", dest="OUTFILE2", help="output file/path")
    parser.add_option("--outfile3", dest="OUTFILE3", help="output file/path")
    
    (options, args) = parser.parse_args()

    out1 = open(options.OUTFILE1,"w")
    out2 = open(options.OUTFILE2,"w")
    out3 = open(options.OUTFILE3,"w")

    out1.write("chr" + "\t" + "pos" + "\t" + "ref" + "\t" + "coverage" + "\t" + "matches" + "\t" + "mismatches" + "\t" + "mismatchRate" + '\n')
    out2.write("chr" + "\t" + "pos" + "\t" + "ref" + "\t" + "coverage" + "\t" + "matches" + "\t" + "mismatches" + "\t" + "mismatchRate" + '\n')
    out3.write("chr" + "\t" + "start" + "\t" + "end" + "\t" + "ref" + "\t" + "coverage" + "\t" + "matches" + "\t" + "mismatches" + "\t" + "mismatchRate" + "\t" + "A" + "\t" + "C" + "\t" + "G" + "\t" + "T" + '\n')
    
    with open('%s' % options.PILEUP) as file:
        next(file)
        for line in file:
            l = str.split(line)
            chr = l[1]
            pos = l[2]
            ref = l[3]
            type = l[4]
            alleles = l[5]
            alleles = list(ref) + alleles.split(",")
            ad = l[8]
            counts = list(map(int,ad.split(",")))
            coverage = sum(counts)
            
            if ((ref == "A") or (ref == "C")):
                if (coverage == 0):
                    matches = 0
                    mismatches1 = 0
                    mismatches2 = 0
                    mmrate1 = 0
                    mmrate2 = 0
                    basecounts={'A':0,'C':0,'G':0,'T':0}
                    
                else:

                    if (type == 'REF'):
                        matches = counts[0]
                        mismatches1 = 0
                        mismatches2 = 0
                        mmrate1 = 0
                        mmrate2 = 0
                        basecounts=basetab(alleles,counts)
                        
                    elif (type == "SNP"):
                        matches = counts[0]
                        mismatches1 = sum(counts[1:])
                        mmrate1 = mismatches1/(matches+mismatches1)
                        mismatches2 = 0
                        mmrate2 = 0
                        basecounts=basetab(alleles,counts)
                    
                    else: # Ignoring INDELs. 
                        continue

                out1.write(chr + "\t" + str(pos) + "\t" + ref + "\t" + str(coverage) + "\t" + str(matches) + "\t" + str(mismatches1) + "\t" + str(mmrate1) + '\n')
                out3.write(chr + "\t" + str(int(pos) -1) + "\t" + str(pos) + "\t" + ref + "\t" + str(coverage) + "\t" + str(matches) + "\t" + str(mismatches1) + "\t" + str(mmrate1) + "\t" + str(basecounts['A']) + "\t" + str(basecounts['C']) + "\t" + str(basecounts['G']) + "\t" + str(basecounts['T']) + '\n')
                
            else: # Reference base must be T or G.
            
                if (coverage == 0):
                    matches = 0
                    mismatches1 = 0
                    mismatches2 = 0
                    mmrate1 = 0
                    mmrate2 = 0
                    basecounts={'A':0,'C':0,'G':0,'T':0}
                    
                else:

                    if (type == 'REF'):
                        matches = counts[0]
                        mismatches1 = 0
                        mismatches2 = 0
                        mmrate1 = 0
                        mmrate2 = 0
                        basecounts=basetab(alleles,counts)
                        
                    elif (type == "SNP"):
                        matches = counts[0]
                        mismatches1 = 0
                        mmrate1 = 0
                        mismatches2 = sum(counts[1:])
                        mmrate2 = mismatches2/(matches+mismatches2)
                        basecounts=basetab(alleles,counts)
                        
                    else: # Ignoring INDELs. 
                        continue                

                out2.write(chr + "\t" + str(pos) + "\t" + ref + "\t" + str(coverage) + "\t" + str(matches) + "\t" + str(mismatches2) + "\t" + str(mmrate2) + '\n')
                out3.write(chr + "\t" + str(int(pos) -1) + "\t" + str(pos) + "\t" + ref + "\t" + str(coverage) + "\t" + str(matches) + "\t" + str(mismatches2) + "\t" + str(mmrate2) + "\t" + str(basecounts['A']) + "\t" + str(basecounts['C']) + "\t" + str(basecounts['G']) + "\t" + str(basecounts['T']) + '\n')
                
    out1.close()
    out2.close()
    out3.close()
