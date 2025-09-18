from optparse import OptionParser
import pandas as pd

# Read a reference table of transcripts and update a query bed file to match.

# Troy Whitfield, 27 Jan., 2023.

if __name__ == '__main__':
    
    parser = OptionParser()
    parser.add_option("--ref", dest="REF", help="Input reference bed file")
    parser.add_option("--query", dest="QUERY", help="Input query bed file")
    parser.add_option("--out", dest="OUT", help="Output bed file")
    
    (options, args) = parser.parse_args()

    bedCols = ['chr','start','end','name','score','strand']
    REF = pd.read_csv(options.REF, sep='\t')
    QUERY = pd.read_csv(options.QUERY, sep='\t', header=None, names=bedCols)
    
    refGene = REF['name']
    refStart = REF['start']
    refEnd = REF['end']
    starts = dict(zip(refGene,refStart))
    ends = dict(zip(refGene,refEnd))
    queryStart = [starts[x] for x in QUERY['name']]
    queryEnd = [ends[x] for x in QUERY['name']]
    QUERY['start'] = queryStart
    QUERY['end'] = queryEnd
    QUERY.to_csv(options.OUT, sep='\t', index=False, header=False)
