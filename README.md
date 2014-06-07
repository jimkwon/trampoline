Trampoline
==========

Trampoline is a set of scripts that generates ready-to-use reference and
annotation data files for general transcriptome-wide researches.


![Logo image](https://raw.githubusercontent.com/hyeshik/trampoline/master/doc/trampoline.png)   
*Modified from [xkcd: Python](http://xkcd.com/353/) by Randall Munroe*


Author
------

Hyeshik Chang <hyeshik@snu.ac.kr>  
Center for RNA Research, Institute for Basic Science, Seoul, South Korea.


Objectives
----------

What is a Trampoline dataset? | What isn't?
----------------------------- | -----------
Single consistent organization of data in standard formats, regardless of genome (or species) | Species-specific layouts, formats, or terms
Ready-to-use bunches for usual cases in transcriptome analysis | Well-organized database that covers all use cases
Simple gene-level analyses | Isoform-level analyses
Rough estimation of expression levels | Perfect resolution of multimapped/unmappable reads
Analysis of gene regulations in known genes | Discovery of novel transcripts or isoforms
Gene structure and terms of multi-cellular eukaryotes | Those of monocellular organism or prokaryotes
Simple text or portable file-based data formats | RDBMS
Requires minimal code for use | Optimized for speed and system resource
Batch processing | Online or real-time processing


Non-redundant RefSeq subset
---------------------------

Trampoline uses a subset of NCBI RefSeq transcripts, called
"non-redundant RefSeq", by filtering out any shorter transcripts
overlapped to longer one in the genome coordinate space. Subsequent
analyses become much easier with the nrRefSeq set. Various structured
databases and bed files provided here enables you to write disposable
transcripts with affordable efforts for explorative analyses. Indeed,
this approach is not suitable when you need isoform-aware processing.


List of components included in the Trampoline dataset
-----------------------------------------------------

* Genome sequences and indices
  - Single uncompressed FASTA format
  - FASTA index (.fai)
  - Chromosome sizes (UCSC Genome Browser format)
  - 2bit sequence (for blat)
* Splice-aware genome sequence indices for:
  - GSNAP/GMAP
  - STAR
* Sequences and indices for highly abundant transcripts
  - Sequence of a rDNA repeat
  - Sequences of 5S, 5.8S rRNA and srpRNA
  - Illumina adapter dimers
  - PhiX genome (which is used for balancing MiSeq runs)
  - Indices of sequences listed above for GSNAP/GMAP and STAR
* Annotations for genomic intervals (in BED format)
  - miRBase hairpins
  - RefSeq transcripts
  - Major RNA classes from Rfam
  - Repeatative sequences from RepeatMasker 
  - tRNAs from gtRNAdb
* Special purpose annotations (in BED format)
  - Transcript end positions
  - nrRefSeq transcripts in BED6 (broken into exons and introns) and BED12 (exon groups).
  - DROSHA and DICER cleavage positions
* Transcript catalogue databases
  - Simple list of nrRefSeq in a text file
  - Transcript sequences in FASTA format
  - Transcript knowledge database in Python 2.x shelve format
  - Transcriptome-genome alignments for super-easy liftOver
    from genome to transcriptome coordinate.


Supported genomes
-----------------

* Human
  - UCSC hg19
  - UCSC hg38
* Mouse *(M. musculus)*
  - UCSC mm10
* Rat *(R. norvegicus)*
  - UCSC rn5
* Pig *(S. scrofa)*
  - UCSC susScr3
* Zebrafish *(D. rerio)*
  - UCSC danRer7
* Fruitfly *(D. melanogaster)*
  - UCSC dm3

