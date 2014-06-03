#!/bin/sh
TMPDIR="$1/rfam-index"
TRNAFULL="$2"
SPECIES="$3"
GENOME2BIT="$4"
OUTPUT="$5"

mkdir -p $TMPDIR

zgrep `echo $SPECIES | sed -e 's, ,.,g'` $2 | \
  awk '{ gsub(/^>/, "", $1); print $1; }' > $TMPDIR/tRNA.sequenceids

zcat $2 | faSomeRecords /dev/stdin $TMPDIR/tRNA.sequenceids $TMPDIR/tRNA.excerpt.fasta

blat -q=rna -minIdentity=90 $GENOME2BIT \
  $TMPDIR/tRNA.excerpt.fasta $TMPDIR/tRNA.psl

pslCDnaFilter -globalNearBest=0.01 -minCover=0.90 $TMPDIR/tRNA.psl \
  $TMPDIR/tRNA.bestHits.psl

pslToBed $TMPDIR/tRNA.bestHits.psl $TMPDIR/tRNA.bed

awk -F'	' 'BEGIN { OFS="\t"; } {
  gsub("^.*trna[0-9]*-", "", $4);
  $4 = sprintf("tRNA|tRNA-%s|tRNA-%s", $4, $4);
  print $0;
}' $TMPDIR/tRNA.bed | cut -f1-6 | bedtools sort | uniq | gzip -c - > $OUTPUT

rm -rf $TMPDIR
