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

TARGETS.extend((
    expand('{genome}/contaminants/contaminants.genomecomp', genome=GENOMES) +
    expand('{genome}/contaminants.star/Genome', genome=GENOMES) +
    expand('{genome}/contaminants.novondx', genome=GENOMES) +
    expand('{genome}/genome/genome.maps/genome.splicesites.iit', genome=GENOMES) +
    expand('{genome}/nrRefSeq-genome.bed.gz', genome=GENOMES) +
    expand('{genome}/annotations.bed.gz', genome=GENOMES) +
    expand('{genome}/nrRefSeq-alignment.db', genome=GENOMES) +
    expand('{genome}/nrRefSeq.fa.fai', genome=GENOMES) +
    expand('{genome}/genome.gff3', genome=GENOMES) +
    expand('{genome}/genome.size', genome=GENOMES) +
    expand('{genome}/genome.star/Genome', genome=GENOMES) +
    expand('{genome}/txEnds.bed.gz', genome=GENOMES) +
    expand('{genome}/nrRefSeq-genome.bed12.gz', genome=GENOMES)
))

rule download_rDNA:
    output: 'downloaded/{genome}/rDNA.fa'
    run:
        url = RDNA_GENBANK_URLS[GENOME2SPECIES[wildcards.genome]]
        shell('wget -O "{output}" "{url}"')

rule rfam_contaminants_seq_names:
    input: 'downloaded/Rfam.fasta.gz'
    output: '{genome}/tmp/Rfam-contaminants-names'
    run:
        species = GENOME2SPECIES[wildcards.genome]
        shell("zgrep '^>.*rRNA.*{species}' {input} | sed -e 's,^>\([^ ]*\).*,\\1,g' > {output}")
        shell("zgrep '^>.*;U3;.*{species}' {input} | sed -e 's,^>\([^ ]*\).*,\\1,g' >> {output}")
        shell("zgrep '^>.*;Metazoa_SRP;.*{species}' {input} | sed -e 's,^>\([^ ]*\).*,\\1,g' >> {output}")

rule make_contaminants_fasta:
    input: rDNA_fasta='downloaded/{genome}/rDNA.fa', \
           rfam_ids='tmp/{genome}/Rfam-contaminants-names', \
           rfam_fasta='downloaded/Rfam.fasta.gz'
    output: '{genome}/contaminants.fa'
    run:
        shell('faSomeRecords {input.rfam_fasta} {input.rfam_ids} {output}.tmp && \
                (cat {input.rDNA_fasta} files/illumina.fa >> {output}.tmp)')
        shell('cat {output}.tmp | awk \'/^>/ {{ print $0; }} \
                    /^[^>]/ {{ gsub(/U/, "T"); print $0; }}\' > {output} && rm -f {output}.tmp')

rule build_contaminants_gsnap_index:
    input: '{genome}/contaminants.fa'
    output: '{genome}/contaminants/contaminants.genomecomp'
    shell: 'gmap_build -T tmp/{wildcards.genome} -D {wildcards.genome} -d contaminants \
                -k 12 -b 12 -q 1 {input}'

rule build_contaminants_star_index:
    input: '{genome}/contaminants.fa'
    output: '{genome}/contaminants.star/Genome'
    threads: 32
    params: outputdir='{genome}/contaminants.star'
    shell: 'STAR --runMode genomeGenerate --genomeDir {params.outputdir} \
                --genomeSAindexNbases 10 --genomeFastaFiles {input} --runThreadN {threads}'

rule build_contaminants_novoalign_index:
    input: '{genome}/contaminants.fa'
    output: '{genome}/contaminants.novondx'
    threads: 32
    shell: 'novoindex -k 10 -s 1 -t 24 {output} {input}'

rule download_genome_sequence:
    output: twobit='{genome}/genome.2bit', fasta='{genome}/genome.fa', \
            fasta_index='{genome}/genome.fa.fai'
    run:
        genome = wildcards.genome
        shell('wget -O {output.twobit} ' + GENOME_2BIT_URL)
        shell('twoBitToFa {output.twobit} {output.fasta}')
        shell('samtools faidx {output.fasta}')

rule build_gsnap_genome_index:
    input: '{genome}/genome.fa'
    output: '{genome}/genome/genome.genomecomp'
    threads: 100 # never run multiple builds
    shell: 'gmap_build -T {wildcards.genome}/downloaded -D {wildcards.genome} \
                -d genome -k 12 -b 12 -q 1 {input}'

rule download_refgene:
    output: 'downloaded/{genome}/refGene.txt.gz'
    run:
        url = REFGENE_URL.format(genome=wildcards.genome)
        shell('wget -O {output} {url}')

rule download_refflat:
    output: 'downloaded/{genome}/refFlat.txt.gz'
    run:
        url = REFFLAT_URL.format(genome=wildcards.genome)
        shell('wget -O {output} {url}')

rule download_reflink:
    output: 'downloaded/{genome}/refLink.txt.gz'
    run:
        url = REFLINK_URL.format(genome=wildcards.genome)
        shell('wget -O {output} {url}')

rule download_knowngene:
    output: 'downloaded/{genome}/knownGene.txt.gz'
    run:
        if wildcards.genome not in NO_KNOWNGENE_AVAILABLE:
            url = KNOWNGENE_URL.format(genome=wildcards.genome)
            shell('wget -O {output} {url}')
        else:
            import gzip
            gzip.open(output[0], 'w')

rule download_refseq_fasta:
    output: 'downloaded/{genome}/refMrna.fa.gz'
    run:
        url = REFMRNA_URL.format(genome=wildcards.genome)
        shell('wget -O {output} {url}')

rule download_refseqali:
    output: 'downloaded/{genome}/refSeqAli.psl'
    run:
        url = REFSEQALI_URL.format(genome=wildcards.genome)
        shell('wget -O - {url} | gzip -d - | sed -n -e "s/^[0-9]*\t//p" > {output}')

rule build_refseqali_db:
    input: alignment='downloaded/{genome}/refSeqAli.psl', \
           nrlist='{genome}/nrRefSeq.list', \
           nrbed='{genome}/nrRefSeq-genome.bed.gz'
    output: '{genome}/nrRefSeq-alignment.db'
    run:
        shell('python tools/build-refSeqAli-database.py {input.nrlist} {input.nrbed} \
                    {input.alignment} {output}')

rule filter_refseq_fasta:
    input: fain='downloaded/{genome}/refMrna.fa.gz', \
           nrlist='{genome}/nrRefSeq.list'
    output: faout='{genome}/nrRefSeq.fa', faidx='{genome}/nrRefSeq.fa.fai'
    run:
        with TemporaryDirectory() as tmpdir:
            refmrnatmp = os.path.join(tmpdir, 'refMrna.fa')
            shell('gzip -cd {input.fain} > {refmrnatmp}')

            shell('faSomeRecords {refmrnatmp} {input.nrlist} {output.faout}')
            shell('samtools faidx {output.faout}')

rule build_gsnap_splice_index:
    input: refgene='downloaded/{genome}/refGene.txt.gz', \
           knowngene='downloaded/{genome}/knownGene.txt.gz', \
           gsnapidx='{genome}/genome/genome.genomecomp'
    output: '{genome}/genome/genome.maps/genome.splicesites.iit'
    shell: '(zcat {input.refgene} | psl_splicesites -s 1; \
             zcat {input.knowngene} | psl_splicesites) | iit_store -o {output}'

rule build_nonredundant_refseq_database:
    input: refflat='downloaded/{genome}/refFlat.txt.gz', \
           reflink='downloaded/{genome}/refLink.txt.gz'
    output: nrdb='{genome}/nrRefSeq.db', nrlist='{genome}/nrRefSeq.list'
    shell: 'zcat {input.refflat} | sort -t "\t" -k3,4 -k5,6n | \
            python tools/build-nonredundant-refseq.py {input.reflink} {output.nrdb} \
              {output.nrlist}'

rule make_nonredundant_refseq_genome_bedanno:
    input: '{genome}/nrRefSeq.db'
    output: '{genome}/nrRefSeq-genome.bed.gz'
    shell: 'python tools/nrrefseq2bed.py {input} | gzip -c - > {output}'

rule make_nonredundant_refseq_genome_bed12anno:
    input: refgene='downloaded/{genome}/refGene.txt.gz', \
           nrreflist='{genome}/nrRefSeq.list'
    output: '{genome}/nrRefSeq-genome.bed12.gz'
    shell: 'python tools/nrrefseq2bed12.py {input.refgene} {input.nrreflist} /dev/stdout | \
                 gzip -c - > {output}'

rule prepare_refseq_catalog:
    input: 'downloaded/{genome}/refGene.txt.gz'
    output: '{genome}/cat.refseq.bed.gz'
    shell: 'python tools/build-refseq-index.py `dirname {output}` {input}'

rule download_miRBase:
    output: 'downloaded/{genome}/mirbase.gff3'
    run:
        URL = MIRBASE_URL.format(species=GENOME2SPECIES_SHORT[wildcards.genome])
        shell("wget -O - {URL} | sed -e 's,^\\([0-9]\\),chr\\1,g' > {output}")

rule prepare_mirbase_catalog:
    input: 'downloaded/{genome}/mirbase.gff3'
    output: '{genome}/cat.mirbase.bed.gz'
    run:
        import urllib, re, gzip

        accession_pat = re.compile('Name=([^;]*)')

        with DeleteOnError(str(output[0]), gzip.open) as out:
            for line in open(input[0]):
                if line.startswith('#'):
                    continue

                fields = line[:-1].split('\t')
                if fields[2] != 'miRNA_primary_transcript':
                    continue

                name = accession_pat.findall(fields[8])[0]
                out.write(('\t'.join([
                    fields[0], str(int(fields[3])-1), fields[4],
                    'miRNA|%s|%s' % (name, name), '.', fields[6]
                ]) + '\n').encode('utf-8'))

rule compile_full_catalog:
    input: '{genome}/cat.mirbase.bed.gz', '{genome}/cat.refseq.bed.gz', \
           '{genome}/cat.rmsk.bed.gz', '{genome}/cat.rfam.bed.gz', \
           '{genome}/cat.trnas.bed.gz'
    output: '{genome}/annotations.bed.gz'
    run:
        inputs = ' '.join(input)
        shell('zcat {inputs} | sort -k1,1 -k2,3n -k4,4 | \
               gzip -c - > {output}')

rule convert_psl_to_gff3:
    input: 'downloaded/{genome}/refSeqAli.psl'
    output: '{genome}/genome.gff3'
    shell: 'tools/blat2gff.pl -match transcript < {input} > {output}'

rule build_star_index:
    input: gff='{genome}/genome.gff3', fasta='{genome}/genome.fa'
    output: '{genome}/genome.star/Genome'
    threads: 32
    params: outputdir='{genome}/genome.star'
    shell: 'STAR --runMode genomeGenerate --genomeDir {params.outputdir} \
                --genomeFastaFiles {input.fasta} --runThreadN {threads} \
                --sjdbGTFtagExonParentTranscript Parent --sjdbGTFfile {input.gff}'

rule make_genome_size:
    input: '{genome}/genome.fa'
    output: '{genome}/genome.size'
    shell: 'faSize -detailed -tab {input} > {output}'

rule make_transcript_ends:
    input: refflat='downloaded/{genome}/refFlat.txt.gz', size='{genome}/genome.size'
    output: '{genome}/txEnds.bed.gz'
    shell: "zcat {input.refflat} | \
            awk -F'	' '{{print $3, $5, $6, $2, 0, $4; }}' | \
            sed -e 's, ,	,g' | \
            bedtools flank -s -l 0 -r 1 -g {input.size} | \
            awk -F'	' '{{ if ($2 < $3) print $0; }}' | \
            bedtools flank -s -l 1 -r 0 -g {input.size} | \
            bedtools sort | gzip -c - > {output}"

# vim: syntax=snakemake
