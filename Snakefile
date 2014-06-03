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

import os
import sys; sys.path.append(os.getcwd())
from itertools import groupby
import csv
from subprocess import PIPE, Popen
from collections import defaultdict

include: 'workflows/snakesupport.py'

REPEATMASKER_IGNORE_CLASSES = 'tRNA snRNA scRNA srpRNA'.split()

GENOMES = ['mm10', 'hg19', 'danRer7', 'dm3']
GENOME2SPECIES = {
    'mm10':     'Mus musculus',
    'hg19':     'Homo sapiens',
    'danRer7':  'Danio rerio',
    'dm3':      'Drosophila melanogaster',
}
GENOME2SPECIES_SHORT = {
    'mm10': 'mmu',
    'hg19': 'hsa',
    'danRer7': 'dre',
    'dm3': 'dme',
}
SUBDIRS = ['downloaded', 'tmp']

#
# Genome-specific data availability
#
NO_KNOWNGENE_AVAILABLE = ['danRer7', 'dm3']
NO_TRNA_IN_UCSC = ['dm3']
SPLIT_REPEAT_MASKER_FILES_IN_UCSC = ['dm3']
CHROMOSOMES = { # for downloading data which are broken into per-chromosome files.
    'dm3': """chr2L chr2LHet chr2R chr2RHet chr3L chr3LHet chr3R chr3RHet
              chr4 chrU chrUextra chrX chrXHet chrYHet chrM""".split(),
}

#
# URLs that are not migrated to new layout yet.
#
GENOME_2BIT_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.2bit'
REFGENE_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/refGene.txt.gz'
REFFLAT_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/refFlat.txt.gz'
REFLINK_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/refLink.txt.gz'
KNOWNGENE_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/knownGene.txt.gz'
MIRBASE_URL = 'ftp://mirbase.org/pub/mirbase/CURRENT/genomes/{species}.gff3'
REFSEQALI_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/database/refSeqAli.txt.gz'
REFMRNA_URL = 'http://hgdownload.soe.ucsc.edu/goldenPath/{genome}/bigZips/refMrna.fa.gz'

GENBANK_EFETCH_URL = (
    'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore'
    '&id={gi}&rettype=fasta&retmode=text')
RDNA_GENBANK_URLS = {
    'Mus musculus': GENBANK_EFETCH_URL.format(gi=38176281),
    'Homo sapiens': GENBANK_EFETCH_URL.format(gi=555853),
    'Danio rerio': GENBANK_EFETCH_URL.format(gi=148614872), # no full rDNA repreat is available
    'Drosophila melanogaster': GENBANK_EFETCH_URL.format(gi=158246),
}

# create directories when missing
for genome in GENOMES:
    for subdir in SUBDIRS:
        dirpath = os.path.join(genome, subdir)
        if not os.path.isdir(dirpath):
            os.makedirs(dirpath)

TARGETS = []

include: 'workflows/legacy.py'
include: 'workflows/RepeatMasker.py'
include: 'workflows/tRNAdb.py'
include: 'workflows/Rfam.py'

rule all:
    input: (lambda _: TARGETS)

# vim: syntax=snakemake
