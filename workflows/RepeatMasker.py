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

REPEATMASKER_FULL_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/rmsk.txt.gz'
REPEATMASKER_SPLIT_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/{chrom}_rmsk.txt.gz'

for genome in GENOMES:
    if genome not in SPLIT_REPEAT_MASKER_FILES_IN_UCSC:
        rule:
            output: 'downloaded/{genome}/rmsk.txt.gz'.format(genome=genome)
            params: genome=genome
            run:
                URL = REPEATMASKER_FULL_URL.format(genome=params.genome)
                shell('wget -O {output} {URL}')
    else:
        rule:
            output: 'downloaded/{genome}/{{chrom}}_rmsk.txt.gz'.format(genome=genome)
            params: genome=genome
            run:
                URL = REPEATMASKER_SPLIT_URL.format(genome=params.genome, chrom=wildcards.chrom)
                shell('wget -O {output} {URL}')

        rule:
            input: expand('downloaded/{genome}/{chrom}_rmsk.txt.gz', genome=[genome], \
                            chrom=CHROMOSOMES[genome])
            output: 'downloaded/{genome}/rmsk.txt.gz'.format(genome=genome)
            shell: 'zcat {input} | gzip -c - > {output}'


rule prepare_repeatmasker_catalog:
    input: 'downloaded/{genome}/rmsk.txt.gz'
    output: '{genome}/cat.rmsk.bed.gz'
    run:
        import gzip

        with DeleteOnError(str(output[0]), gzip.open) as out:
            for line in gzip.open(str(input[0])):
                fields = line[:-1].split(b'\t')
                if fields[12] in REPEATMASKER_IGNORE_CLASSES:
                    continue

                out.write(b'\t'.join([
                    fields[5], fields[6], fields[7],
                    fields[11] + b'|' + fields[10] + b'|' + fields[12] + b'/' + fields[10],
                    fields[1], fields[9]]) + b'\n')

# vim: syntax=snakemake
