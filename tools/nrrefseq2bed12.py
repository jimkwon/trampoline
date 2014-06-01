import pandas as pd
import sys

refgene_file = sys.argv[1]
nrrefseq_file = sys.argv[2]
output_file = sys.argv[3]

def adopt_comma_split_ints(s):
    return map(int, s.rstrip(',').split(','))

REFGENE_COLUMNS = [
    'bin', 'name', 'chrom', 'strand', 'txStart', 'txEnd',
    'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'score',
    'name2', 'cdsStartStat', 'cdsEndStat', 'exonFrames',
]

refgene = pd.read_table(refgene_file, names=REFGENE_COLUMNS, compression='gzip')
nrrefseq = pd.read_table(nrrefseq_file, names=('accession',))

merged = pd.merge(nrrefseq, refgene, how='inner', left_on='accession', right_on='name')
merged['thickStart'] = merged['txStart']
merged['thickEnd'] = merged['txEnd']
merged['itemRgb'] = '0,0,0'
merged['blockSizes'] = [','.join(str(end - start)
                         for start, end in zip(adopt_comma_split_ints(row['exonStarts']),
                                               adopt_comma_split_ints(row['exonEnds'])))
                        for rowi, row in merged.iterrows()]
merged['blockStarts'] = [row['exonStarts'].rstrip(',') for rowi, row in merged.iterrows()]

bed12 = merged[['chrom', 'txStart', 'txEnd', 'name', 'score', 'strand', 'thickStart', 'thickEnd', 'itemRgb', 'exonCount', 'blockSizes', 'blockStarts']]
bed12.sort(['chrom', 'txStart']).to_csv(output_file, sep='\t', header=False, index=False)
