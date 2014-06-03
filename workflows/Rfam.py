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

RFAM_FASTA_URL = 'ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.fasta.gz'
RFAM_FULL_URL = 'ftp://ftp.sanger.ac.uk/pub/databases/Rfam/CURRENT/Rfam.full.gz'

rule download_rfam_fasta:
    output: 'downloaded/Rfam.fasta.gz'
    shell: 'wget -O {output} {RFAM_FASTA_URL}'

rule download_rfam_full:
    output: 'downloaded/Rfam.full.gz'
    shell: 'wget -O {output} {RFAM_FULL_URL}'

rule prepare_rfam_catalog:
    input: rfamfasta='downloaded/Rfam.fasta.gz', \
           rfamfull='downloaded/Rfam.full.gz', \
           twobit='{genome}/genome.2bit'
    output: '{genome}/cat.rfam.bed.gz'
    run:
        species = GENOME2SPECIES[wildcards.genome]
        with TemporaryDirectory() as tmpdir:
            shell('sh tools/build-rfam-index.sh {tmpdir} {input.rfamfasta} ' \
                  '{input.rfamfull} "{species}" {input.twobit} {output}')

# vim: syntax=snakemake
