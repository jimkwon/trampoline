#
# Copyright (c) 2011-2014 Hyeshik Chang
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all 
# copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE 
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE 
# SOFTWARE.
#

TRNADB_FULL_URL = 'http://gtrnadb.ucsc.edu/download/tRNAs/GtRNAdb-all-tRNAs.fa.gz'
TRNADB_UCSC_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/tRNAs.txt.gz'

rule download_gtRNAdb_sequences:
    output: 'downloaded/gtRNAdb-all.fa.gz'
    run:
        URL = TRNADB_FULL_URL
        shell('wget -O {output} {URL}')

rule download_tRNAdb_from_UCSC:
    output: 'downloaded/{genome}/tRNAs.txt.gz'
    run:
        URL = TRNADB_UCSC_URL.format(genome=wildcards.genome)
        shell('wget -O {output} {URL}')


for genome in GENOMES:
    if genome not in NO_TRNA_IN_UCSC:
        rule: # prepare_trna_catalog
            input: 'downloaded/{genome}/tRNAs.txt.gz'.format(genome=genome)
            output: '{genome}/cat.trnas.bed.gz'.format(genome=genome)
            run:
                import re, gzip

                with DeleteOnError(str(output[0]), gzip.open) as out:
                    for line in gzip.open(str(input[0])):
                        fields = line[:-1].split(b'\t')
                        out.write(b'\t'.join([
                            fields[1], fields[2], fields[3],
                            b'tRNA|tRNA-' + fields[7] + fields[8] + b'|tRNA-' +
                            fields[7] + fields[8],
                            fields[5], fields[6]]) + b'\n')
    else:
        rule:
            input: trna='downloaded/gtRNAdb-all.fa.gz', \
                   twobit='{genome}/genome.2bit'.format(genome=genome)
            output: '{genome}/cat.trnas.bed.gz'.format(genome=genome)
            params: genome=genome
            run:
                species = GENOME2SPECIES[params.genome]
                with TemporaryDirectory() as tmpdir:
                    shell('sh tools/build-trna-index.sh {tmpdir} {input.trna} ' \
                          '"{species}" {input.twobit} {output}')

# vim: syntax=snakemake
